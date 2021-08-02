
R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
Copyright (C) 2021 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> setwd("~/2021-REU/SNV")
> #BiocManager::install("VariantAnnotation")
> library(VariantAnnotation)
> library(GenomicRanges)
> library(tidyverse)
> library(janitor)
> library(GenomicFeatures)
> load("~/2021-REU/CNV Analysis/partials/sample_data")
> sample_data_all <-
+   read.csv("~/2021-REU/samplelist.csv") %>% clean_names()%>% mutate(sample_name = str_replace_all(sample_name, " - HiSeq", "")) %>% mutate(type = grepl("FP", tissue)) %>% mutate(sample_name = str_split(sample_name, " ",n = 1))
> 
> #ref<- getChromInfoFromNCBI("GCF_015237465.1")$RefSeqAccn
> genes <- genes(loadDb("~/2021-REU/CNV Analysis/rCheMyd1.sqlite"))
> #genes_by_chr <- splitAsList(genes, seqnames(genes))
> #rm(genes)
> fl <- "inputs/Somatypus_Indels_final.vcf.gz"
> 
> # need to load one chromosome at a time or its too big! + currently subsetting to exons
> # only need the NR and NV data from GENO, nothing from INFO
> 
> convertToMatrix <- function(x){
+   rows <- rownames(x)
+   cols <- colnames(x)
+   data <- as.numeric(x)
+   dim(data) <- c(length(rows),length(cols))
+   dimnames(data) <- list(rows, cols)
+   return(data)
+ }
> 
> ### PLOTTING FUNCTIONS
> getX <- function(ind, tumor=TRUE){ # getX for VAF
+   variants <- names(vars.vaf[tumour.only.idx == tumor, ind])
+   x <- ranges[names(ranges) %in% variants]
+   x <- x[match(variants, names(x)),]
+   return(start(x))
+ }
> 
> plotVAF <- function(ind,ylab=NA){ # plot the vaf
+   # plot the normals
+   plot(getX(ind, tumor=FALSE), vars.vaf[!tumour.only.idx, ind],
+        ylab=ylab, xlab=names[ind], pch=20, col="gray72", cex=0.6,xlim=xlim)
+   # then the tumor in red
+   points(getX(ind, tumor=TRUE), vars.vaf[tumour.only.idx, ind[1]], pch=18, col="red", cex=0.7)
+ }
> 
> plotCoverage <- function(ind,ylab=NA){
+   # Coverage (log10): all variant
+   x <- ranges
+   y <- log10(vars.nr[,ind])
+   plot(start(x[match(names(y), names(x)),]), y, ylab=ylab,
+        xlab=names[ind], pch=20, col="gray72", cex=0.7,ylim=c(0,2),xlim=xlim)
+ }
> 
> plotFull <- function(toPlot,name){
+   i <- seq_along(samples)[toPlot]
+   # Create PNG file
+   png(paste0(length(vars.vaf[tumour.only.idx,]), "_", names(gene), "_", name, ".png"), 2400, 1000)
+   par(mfrow=c(2,length(samples[toPlot])))
+   plotVAF(i[1],ylab="VAF")
+   mtext(paste("VAF for gene", names(gene), "in sample"), outer = TRUE, side=3,line=-2)
+   for (j in i[2:length(i)]) {
+     plotVAF(j)
+   }
+ 
+   plotCoverage(i[1],ylab="log10(Coverage)")
+ 
+   for (j in i[2:length(i)]) {
+     plotCoverage(j)
+     # Close PNG
+   }
+   mtext(paste("Coverage for gene", names(gene), "in sample"),outer = TRUE, side=3,line=-55)
+   dev.off()
+ }
> 
> start_gene <- 304
> display_percent <- 100
> for(g in seq_along(genes)[start_gene:length(genes)]){
+   gene <- genes[g]
+   curr_percent <-round(g/length(genes) * 100,digits=2)
+   if(display_percent != curr_percent){
+   message(names(gene), " - ",g, "/",length(genes), "(", curr_percent, "%)")
+     display_percent <- curr_percent
+   }
+   gc()
+   setwd("~/2021-REU/SNV")
+   tryCatch({
+     svp <- ScanVcfParam(geno=c("NR","NV"), info=NA, which=gene)
+     vcf <- readVcf(fl,as.character(seqnames(gene)), svp)
+ 
+     vars.nv <- convertToMatrix(geno(vcf)$NV)
+     if(length(vars.nv) > 1){
+       vars.nr <- convertToMatrix(geno(vcf)$NR)
+ 
+       vars.vaf <- vars.nv/vars.nr
+ 
+       ### Get samples
+       samples <- colnames(vars.nv)
+       names <- str_remove_all(samples, "bwa_|_CheMyd|_sorted")
+       tumours <- grepl(paste((sample_data_all %>% filter(type))$sample_name, sample_data$sample_name, sep="|",  collapse="|"), samples,ignore.case = TRUE) & !grepl("HiSeq", samples, ignore.case=TRUE)
+       hosts <- grepl(paste((sample_data_all %>% filter(!type))$sample_name, sample_data$normal, sep="|", collapse="|"), samples,ignore.case=TRUE) & !grepl("HiSeq", samples, ignore.case=TRUE)
+       samples[!(tumours | hosts)] # make sure we've labeled everything except HiSeq
+ 
+ 
+       ### BASIC MANIPULATION
+       MIN.VAF <- 0.25
+       not.hosts <- apply(vars.vaf[,hosts], 1, function(vaf) {
+         all(vaf < MIN.VAF)
+       })
+ 
+       # Second, take the variants that are found with at least 3 supporting reads in at least one tumour
+       MIN.READS <- 3
+       one.tumour <- apply(vars.nv[,tumours], 1, function(nv) {
+         any(nv >= MIN.READS)
+       })
+ 
+       # Extract tumour-only variants and host variants
+       tumour.only.idx <- not.hosts & one.tumour
+       tumour.only.idx[is.na(tumour.only.idx)] <- FALSE # when in doubt, mark false (host)
+ 
+       tumour.only.snvs <- rowRanges(vcf)[tumour.only.idx,]
+       host.snvs <- rowRanges(vcf)[!tumour.only.idx,]
+       #rm(vcf)
+ 
+       # 3) Take histogram showing number of samples the T-only variants are shared between
+ 
+       # Get number of samples in which each variant is present (≥3 reads)
+       #num.tumours <- apply(vars.nv[tumour.only.idx, tumours], 1, function(nv) {
+       #  sum(nv >= MIN.READS)
+       #})
+       #hist(num.tumours, breaks=100, col="cornflowerblue", border="white", main="Number of tumours containing each T-only variant (≥3 reads)", xlab=NULL)
+ 
+ 
+       setwd("outputs/BasicManipulation/")
+       xlim <- c(start(range(gene)), end(range(gene)))
+       ranges <-ranges(rowRanges(vcf))
+ 
+       plotFull(tumours, "tumor")
+       plotFull(hosts,"host")
+       rm(vars.nr)
+       rm(vars.nv)
+       rm(vars.vaf)
+     }
+   },
+   error=function(cond) {
+     print(paste("ERROR: Problem with gene", names(gene), "number", g))
+     print(cond)
+     list <- dev.list()
+     tryCatch(dev.off(),error=function(w){print("Nothing to close.")})
+   })
+ }
[1] "ERROR: Problem with gene ADRB1 number 307"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene ADRM1 number 310"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AEBP1 number 314"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AEBP2 number 315"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AFTPH number 327"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AGAP2 number 330"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AGBL3 number 334"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AGL number 341"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AGPAT4 number 348"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AGT number 354"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AHCTF1 number 360"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AHNAK number 366"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AIF1L number 374"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AK1 number 383"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AK4 number 386"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AKAP11 number 394"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AKR1A1 number 410"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AKR1D1 number 411"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AKT1S1 number 414"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ALG1 number 438"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ALG11 number 439"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ALG12 number 440"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ALG6 number 445"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ALG8 number 446"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ALG9 number 447"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ALK number 448"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ALKAL1 number 449"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ALKBH8 number 458"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ALLC number 459"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ALPK3 number 464"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ALS2CL number 467"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ALX4 number 470"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ALYREF number 471"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AMDHD1 number 476"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AMER2 number 478"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AMH number 481"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AMHR2 number 482"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AMIGO3 number 485"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene AMPD2 number 494"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ANAPC10 number 501"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ANGPT1 number 512"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ANKDD1A number 526"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ANKRA2 number 538"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ANKRD12 number 542"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ANKRD13A number 543"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ANKRD13C number 545"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ANKRD16 number 547"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ANKRD2 number 549"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ANKRD33 number 558"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene ANKRD34C number 562"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene ANKRD45 number 569"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ANKRD6 number 577"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ANKS3 number 584"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ANO6 number 596"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ANXA11 number 609"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ANXA13 number 610"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene AP1M2 number 627"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AP3B1 number 634"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AP3B2 number 635"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AP3M1 number 637"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AP5B1 number 644"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene AP5S1 number 646"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene APAF1 number 648"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene APBB1 number 652"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene APBB3 number 654"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene APC2 number 656"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene APCDD1 number 657"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene APELA number 660"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene API5 number 665"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene APLF number 666"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene APLP1 number 669"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene APOA1 number 671"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene APOA5 number 674"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene APOH number 683"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene APPBP2 number 688"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene APPL1 number 689"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARAP3 number 705"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AREL1 number 709"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARF1 number 710"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARFGAP3 number 717"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARFGEF2 number 719"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARFRP1 number 723"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARG2 number 725"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARHGAP29 number 743"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARHGAP31 number 745"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARHGAP35 number 748"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARHGAP44 number 753"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARHGEF1 number 762"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARHGEF39 number 779"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARHGEF6 number 782"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARHGEF7 number 783"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARID3A number 788"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene ARL4D number 810"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene ARMC10 number 823"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARMT1 number 836"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARPC4 number 845"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARPP19 number 849"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARRDC2 number 854"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARRDC4 number 856"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ARSJ number 862"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ASB1 number 875"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ASB16 number 882"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ASCL3 number 897"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ASIC5 number 906"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ASL number 907"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ASMTL number 908"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATAD1 number 929"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATF6B number 943"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618392.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene ATG101 number 947"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATG16L2 number 952"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATG4B number 957"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATG4C number 958"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATG5 number 960"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATG9A number 962"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATL3 number 966"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATOH8 number 972"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATOX1 number 973"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP10B number 975"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP10D number 976"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP11B number 978"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP11C number 979"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP1B1 number 987"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP2A2 number 992"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP2C1 number 997"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP4A number 999"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP5F1E number 1005"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP5MK number 1013"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP5PF number 1016"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP6AP2 number 1020"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP6V0A2 number 1022"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP6V0A4 number 1023"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP6V1B2 number 1031"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP7B number 1042"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP8B3 number 1047"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP8B4 number 1048"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATP9B number 1050"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATPSCKMT number 1053"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATXN10 number 1060"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ATXN1L number 1061"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene ATXN2L number 1063"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AVL9 number 1077"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AVPI1 number 1078"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AXIN2 number 1083"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene AZI2 number 1085"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene B3GALNT1 number 1089"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene B3GALT2 number 1092"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene B4GALT2 number 1112"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene B4GALT3 number 1113"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene B4GAT1 number 1118"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BAALC number 1121"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BAAT number 1122"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BACE2 number 1126"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BACH1 number 1127"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BAG2 number 1131"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BAG5 number 1134"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BANF1 number 1142"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BASP1 number 1151"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BATF2 number 1153"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene BAZ1B number 1157"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BAZ2A number 1158"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BCHE number 1183"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BCKDK number 1186"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BCL11A number 1188"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BCL2A1 number 1191"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene BCL2L10 number 1193"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene BCL2L13 number 1196"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BCL2L14 number 1197"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BCL2L15 number 1198"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene BCL7B number 1202"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene BCLAF1 number 1205"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BCR number 1211"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BDNF number 1215"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BDP1 number 1216"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BEND3 number 1220"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BEND6 number 1223"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BGN number 1234"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BHLHE22 number 1237"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene BICD2 number 1244"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BICRAL number 1248"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BIN3 number 1252"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BMP1 number 1271"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BMPR1A number 1282"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BMX number 1287"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BNC1 number 1288"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BNIP3 number 1292"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BOD1L1 number 1298"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BOLA3 number 1302"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BORCS6 number 1307"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene BORCS7 number 1308"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BPHL number 1310"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BPNT1 number 1314"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BRAF number 1317"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BRAT1 number 1318"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene BRCA1 number 1319"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BRD8 number 1328"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BRIP1 number 1338"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BRMS1L number 1342"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BTBD1 number 1360"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BTBD3 number 1368"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BTK number 1379"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene BZW1 number 1390"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C10H16orf89 number 1399"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C11H2orf69 number 1405"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C12H17orf67 number 1409"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene C14H16orf46 number 1418"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C18H1orf167 number 1437"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C1H11orf54 number 1448"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C1H12orf4 number 1453"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C1H12orf40 number 1454"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C1H12orf50 number 1456"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C1H12orf57 number 1457"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene C1H21orf62 number 1462"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C1H3orf52 number 1467"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C1H3orf85 number 1468"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C1HXorf38 number 1469"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C1S number 1488"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene C24H6orf89 number 1504"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C2CD2 number 1513"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C2H18orf21 number 1521"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C2H6orf201 number 1526"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C2H6orf62 number 1528"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C3H1orf131 number 1540"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C3H6orf120 number 1545"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene C4H20orf194 number 1549"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene C4H4orf45 number 1555"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C4H4orf47 number 1557"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C4H4orf48 number 1558"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C5 number 1561"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C5H11orf16 number 1563"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C5H11orf24 number 1564"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C5H14orf39 number 1572"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C5H15orf62 number 1573"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene C6H9orf135 number 1581"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C7H10orf143 number 1588"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C7H10orf90 number 1593"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C7H11orf98 number 1599"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene C7H3orf14 number 1600"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C7H3orf18 number 1601"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C8B number 1605"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene C8H1orf112 number 1607"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C8H5orf24 number 1615"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene C9HXorf65 number 1624"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene C9orf72 number 1625"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CA14 number 1628"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CAB39 number 1636"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CABYR number 1644"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CACNA1H number 1654"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CALCOCO1 number 1686"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CALCOCO2 number 1687"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CALD1 number 1690"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CALHM6 number 1694"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene CALML6 number 1700"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CALU number 1704"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CAMK1 number 1706"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CAMK2B number 1710"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CAPN10 number 1733"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CAPN11 number 1734"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CAPNS1 number 1744"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CASD1 number 1770"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CASKIN1 number 1772"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CASKIN2 number 1773"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CASP8AP2 number 1779"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CATSPER4 number 1792"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CATSPERG number 1796"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CAV3 number 1799"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CAVIN1 number 1800"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CAVIN2 number 1801"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CBL number 1806"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CC2D2B number 1829"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC105 number 1836"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC110 number 1839"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene CCDC113 number 1841"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC12 number 1845"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC120 number 1846"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC127 number 1851"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC134 number 1854"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC136 number 1855"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC137 number 1856"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC14 number 1858"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC142 number 1860"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC152 number 1866"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC159 number 1870"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC166 number 1872"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene CCDC181 number 1885"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC186 number 1889"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC190 number 1892"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC28B number 1903"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC33 number 1907"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene CCDC34 number 1908"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC40 number 1911"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC51 number 1916"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC59 number 1918"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC6 number 1919"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC61 number 1921"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC62 number 1922"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC66 number 1925"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC73 number 1931"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC77 number 1933"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCDC85B number 1940"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene CCDC86 number 1942"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCN4 number 1971"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCN5 number 1972"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCNA1 number 1974"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCNB2 number 1978"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCND2 number 1982"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCND3 number 1983"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCNJL number 1993"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCNL1 number 1995"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCT5 number 2014"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CCT8 number 2017"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CD101 number 2019"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CD164L2 number 2023"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CD226 number 2029"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CD244 number 2030"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CD248 number 2032"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene CD4 number 2044"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CD47 number 2048"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CD58 number 2051"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CD68 number 2056"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618375.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene CD82 number 2062"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CD8A number 2064"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CD96 number 2068"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CD99L2 number 2069"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDC123 number 2073"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDC20 number 2077"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDC20B number 2078"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDC26 number 2083"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDC37L1 number 2087"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDC42EP5 number 2097"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDC7 number 2103"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDH10 number 2114"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDHR5 number 2137"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDK15 number 2146"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDK5R1 number 2158"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene CDK5RAP1 number 2160"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDKN2AIPNL number 2174"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDKN3 number 2178"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDNF number 2179"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDRT4 number 2186"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDS1 number 2187"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDS2 number 2188"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDV3 number 2190"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDX1 number 2191"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CDX4 number 2193"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CEBPA number 2197"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene CEBPZ number 2202"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CELF3 number 2208"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CELF5 number 2210"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CELF6 number 2211"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CENPI number 2222"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CENPL number 2225"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CENPT number 2230"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CEP120 number 2237"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CEP131 number 2240"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CEP152 number 2242"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CEP164 number 2244"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CEP290 number 2250"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CEP57L1 number 2258"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CEP63 number 2259"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CEP83 number 2265"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CERK number 2274"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CERS3 number 2278"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CFAP161 number 2286"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CFAP20 number 2287"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CFAP20DC number 2289"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CFAP299 number 2293"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CFAP36 number 2295"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CFAP65 number 2307"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CFAP91 number 2313"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CFDP1 number 2320"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CFL2 number 2324"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CGA number 2327"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CGREF1 number 2332"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHADL number 2338"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHCHD10 number 2344"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHCHD7 number 2350"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHD1 number 2351"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHD3 number 2354"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHD5 number 2356"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHD8 number 2359"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHIC2 number 2368"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHMP4B number 2378"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHN1 number 2383"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHODL number 2385"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHORDC1 number 2386"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHPF number 2388"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHRDL2 number 2394"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHRM5 number 2400"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene CHRNA3 number 2404"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHRNA4 number 2405"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHRNA9 number 2408"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene CHRNB1 number 2409"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618375.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene CHRNB3 number 2411"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHRNB4 number 2412"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHRNE number 2414"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHST13 number 2420"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHTF18 number 2431"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CHUK number 2433"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CIAO1 number 2435"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CIB3 number 2443"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CIBAR2 number 2445"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CIDEC number 2449"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CIP2A number 2453"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CIR1 number 2455"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CKAP2 number 2465"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CKAP4 number 2467"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CKM number 2470"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CLDN1 number 2485"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CLDN10 number 2486"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CLDN12 number 2488"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CLDN14 number 2489"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene CLDN18 number 2492"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CLDN25 number 2497"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene CLEC3A number 2504"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene CLIC6 number 2514"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CLK4 number 2522"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CLN3 number 2524"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CLN6 number 2526"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CLPTM1 number 2533"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CLSTN1 number 2540"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CLUL1 number 2549"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CMC1 number 2554"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CMC2 number 2555"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CMTM5 number 2563"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CNGA4 number 2578"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CNIH3 number 2583"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CNKSR1 number 2585"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CNMD number 2588"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CNNM3 number 2594"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CNOT3 number 2600"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CNOT4 number 2601"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CNOT6 number 2602"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CNOT9 number 2606"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CNPPD1 number 2607"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CNTRL number 2627"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COBLL1 number 2631"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COG1 number 2633"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COG2 number 2634"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COG7 number 2639"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COL10A1 number 2642"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COL12A1 number 2645"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COL16A1 number 2649"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COL28A1 number 2663"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COL2A1 number 2664"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COL3A1 number 2665"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COL5A2 number 2672"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COL6A1 number 2674"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COL6A6 number 2677"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COL9A2 number 2682"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COLGALT2 number 2689"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COMMD1 number 2691"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COMMD3 number 2694"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene COMT number 2700"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COMTD1 number 2701"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COPZ1 number 2719"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COPZ2 number 2720"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene COQ2 number 2723"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CORO1A number 2733"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CPB2 number 2743"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CPE number 2745"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CPNE9 number 2767"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CPOX number 2769"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CPS1 number 2772"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CPSF1 number 2773"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CPT2 number 2782"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CRABP1 number 2787"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CRABP2 number 2788"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CRB3 number 2798"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CREB1 number 2800"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CREB3L2 number 2803"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CREB3L4 number 2805"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CREM number 2814"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CRH number 2815"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene CRISPLD2 number 2823"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CRK number 2824"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CRLF1 number 2826"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CRLF2 number 2827"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CRYAA number 2841"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CRYZL1 number 2856"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CSAD number 2858"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CSDE1 number 2860"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CSF1 number 2862"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CSF3R number 2865"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CSNK1D number 2874"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CSNK2B number 2880"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CSRNP2 number 2885"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CSTF3 number 2894"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CTCF number 2899"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CTCFL number 2900"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CTDSPL number 2904"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CTDSPL2 number 2905"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CTH number 2907"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CTNNA1 number 2910"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CTPS1 number 2920"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CTR9 number 2922"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CTSA number 2924"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CTSK number 2931"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CTSL number 2932"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CTSO number 2933"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CTTNBP2 number 2937"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CTXN2 number 2941"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene CWC15 number 2958"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CWH43 number 2964"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CXCL12 number 2967"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CYGB number 2978"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CYREN number 2983"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CYRIA number 2984"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CYRIB number 2985"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene CZIB number 2994"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DAB2 number 2999"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DAO number 3010"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DAP3 number 3011"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DAZAP1 number 3021"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DAZAP2 number 3022"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DBN1 number 3028"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DBNDD1 number 3029"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DBNL number 3031"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DCAF12 number 3039"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DCP1A number 3066"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DCP1B number 3067"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DCT number 3072"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DCTN5 number 3078"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DDIT4L number 3099"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene DDO number 3100"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DDR1 number 3102"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DDX17 number 3108"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DDX27 number 3114"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene DDX3X number 3118"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DDX42 number 3121"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DDX54 number 3129"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DDX6 number 3134"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DECR1 number 3136"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DEDD number 3138"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DENND4B number 3157"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DENND6A number 3161"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DEPP1 number 3169"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene DESI1 number 3177"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DESI2 number 3178"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DET1 number 3179"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DFFA number 3182"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DGCR2 number 3186"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DGKA number 3188"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DGKD number 3190"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DHFR number 3204"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DHRS1 number 3208"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene DHRS13 number 3211"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene DHRS7C number 3215"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DHX57 number 3230"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DHX9 number 3233"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DIAPH1 number 3235"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DIPK2A number 3250"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DIRAS1 number 3252"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DIS3L number 3256"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DLEC1 number 3271"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DLK2 number 3283"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DLX3 number 3290"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DLX5 number 3292"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene DLX6 number 3293"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DMAC1 number 3294"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DMC1 number 3300"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DMGDH number 3302"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DMTN number 3312"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DNA2 number 3316"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DNAAF3 number 3319"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DNAAF6 number 3322"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DNAH6 number 3332"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DNAI2 number 3337"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DNAJA2 number 3341"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DNAJB11 number 3345"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DNAJB12 number 3346"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DNAJB5 number 3350"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DNAJC10 number 3355"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DNAJC16 number 3361"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DNAJC25 number 3369"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DNAJC5G number 3377"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DNASE2B number 3387"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DNM1L number 3392"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DNMT3B number 3397"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DOK4 number 3419"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DOLPP1 number 3424"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DPF3 number 3434"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DPP8 number 3449"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DPP9 number 3450"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DPY19L1 number 3452"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DPY19L4 number 3454"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DRAM2 number 3465"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DRD3 number 3473"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DRGX number 3478"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DRP2 number 3480"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DSEL number 3484"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DSN1 number 3485"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DTD2 number 3490"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DTHD1 number 3491"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DUS3L number 3507"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DUSP10 number 3510"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DUSP14 number 3513"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DUSP15 number 3514"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DUSP18 number 3516"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene DUSP19 number 3517"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DUSP5 number 3524"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DUSP8 number 3527"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DVL1 number 3529"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DVL3 number 3530"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DYNC1H1 number 3532"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DYNC2I1 number 3538"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DYNLT2 number 3546"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DYRK4 number 3554"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DZANK1 number 3557"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene DZIP1L number 3559"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene E2F2 number 3561"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EAF1 number 3569"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EAF2 number 3570"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene EARS2 number 3571"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EBAG9 number 3572"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EBF3 number 3575"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ECE1 number 3581"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ECEL1 number 3583"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ECI2 number 3590"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ECM2 number 3592"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ECSCR number 3595"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EDA number 3599"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EDARADD number 3602"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EDN1 number 3609"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EEF1AKMT1 number 3619"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EEF1E1 number 3626"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EEF2 number 3628"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EEF2K number 3629"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EFCAB7 number 3642"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EFCAB8 number 3643"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EFEMP1 number 3646"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EFNA2 number 3655"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EGF number 3666"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene EGLN3 number 3672"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EHD3 number 3681"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EI24 number 3686"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EIF1AX number 3689"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EIF2AK4 number 3695"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EIF2B5 number 3700"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EIF2D number 3701"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EIF3G number 3710"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EIF3H number 3711"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EIF4A1 number 3717"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618375.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene EIF4A3 number 3719"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EIF4E number 3721"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EIF4E3 number 3724"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EIF4EBP1 number 3725"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EIF4EBP3 number 3727"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ELAC1 number 3738"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ELAC2 number 3739"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ELL2 number 3756"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ELOA number 3765"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene ELOF1 number 3768"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene ELP1 number 3776"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ELP5 number 3780"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618351.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene EMC3 number 3785"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EMCN number 3791"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EML1 number 3799"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EML3 number 3801"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene ENC1 number 3813"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ENDOU number 3816"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ENDOV number 3817"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ENPP4 number 3833"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ENPP5 number 3834"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ENPP6 number 3835"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EPB41 number 3854"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EPC1 number 3860"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EPHA1 number 3866"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EPHA2 number 3868"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EPN3 number 3885"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ERBB3 number 3901"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ERCC5 number 3910"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ERH number 3922"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ERN2 number 3939"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ERP27 number 3942"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene ESCO1 number 3947"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ESPL1 number 3952"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ESPN number 3953"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ESS2 number 3961"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ETFA number 3968"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ETS1 number 3975"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ETV3 number 3978"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ETV7 number 3983"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EVPL number 3994"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EVX2 number 3996"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene EXOC3L1 number 4007"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EXOC3L2 number 4008"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EXOSC10 number 4018"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EXT2 number 4028"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene EZR number 4039"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene F8 number 4052"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAAP100 number 4056"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FABP3 number 4061"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FADS6 number 4066"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAM155B number 4110"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAM162A number 4117"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAM166C number 4123"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAM167A number 4124"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAM174A number 4134"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAM183A number 4143"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAM189B number 4150"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAM199X number 4153"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAM214A number 4161"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FAM217A number 4164"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAM221A number 4169"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAM222B number 4172"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAM229A number 4176"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FAM234A number 4177"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAM32A number 4184"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FAM43A number 4189"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FAM72A number 4197"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAM76B number 4199"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAM83C number 4206"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAM83E number 4207"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAM98A number 4215"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAM98B number 4216"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FANCG number 4226"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FANCL number 4228"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FANK1 number 4230"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FASTK number 4241"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FASTKD3 number 4244"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FAT2 number 4247"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FBLIM1 number 4256"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FBLN2 number 4258"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FBLN5 number 4259"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FBXL18 number 4273"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FBXO30 number 4296"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FBXO41 number 4305"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FBXO5 number 4312"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FBXO7 number 4313"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FCHO1 number 4326"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FDX1 number 4334"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FDXR number 4337"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FEZ2 number 4354"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FGB number 4359"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FGD2 number 4361"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FGD5 number 4364"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FGF1 number 4366"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FGF11 number 4368"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618375.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FGF20 number 4377"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FGF22 number 4379"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FGF23 number 4380"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FGF3 number 4381"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FGF4 number 4382"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FGF5 number 4383"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FGF7 number 4385"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FGFBP3 number 4390"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FGFR3 number 4394"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FGGY number 4398"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FGL1 number 4399"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FHL3 number 4408"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FIGN number 4417"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FILIP1L number 4421"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FIZ1 number 4426"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FKBP14 number 4430"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FKBP6 number 4438"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FKBP9 number 4441"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FKBPL number 4442"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618392.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FLACC1 number 4444"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FLI1 number 4447"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FLRT2 number 4454"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FLVCR1 number 4459"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FMC1 number 4461"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FMN2 number 4463"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FMO2 number 4467"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FMR1NB number 4471"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FNBP1 number 4474"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FNBP4 number 4476"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FNTA number 4488"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FOXB2 number 4498"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FOXC1 number 4499"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FOXD4L3 number 4504"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene FOXE3 number 4506"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene FOXJ2 number 4515"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FOXK1 number 4516"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FOXL3 number 4520"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FOXN1 number 4522"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FOXO4 number 4528"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FOXP4 number 4533"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FOXRED1 number 4536"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FPGS number 4539"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FRAT2 number 4542"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FREM3 number 4545"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FRMPD1 number 4556"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FSD1 number 4569"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FSHB number 4572"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FSTL3 number 4578"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FTCD number 4581"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FTH1 number 4583"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FUBP3 number 4588"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FUCA2 number 4590"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FUT11 number 4596"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FUT7 number 4598"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FZD1 number 4612"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FZD10 number 4613"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FZD3 number 4615"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FZD7 number 4619"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene FZD8 number 4620"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene FZR1 number 4622"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene G3BP2 number 4626"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene G6PC3 number 4628"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GAA number 4630"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GABPA number 4637"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GABPB1 number 4638"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GABRA2 number 4641"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GABRA3 number 4642"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GABRB3 number 4648"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GABRQ number 4652"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene GABRR3 number 4654"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GADD45G number 4659"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene GALNT6 number 4688"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GALR2 number 4694"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GAPDH number 4702"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GART number 4710"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GAS2L3 number 4715"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GAS8 number 4718"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GASK1B number 4720"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GATA3 number 4723"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GATA6 number 4726"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GATAD2A number 4728"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GATB number 4730"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GATD1 number 4732"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GBX2 number 4740"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene GDF2 number 4774"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene GDF7 number 4777"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GDI2 number 4779"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GDPD2 number 4782"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GDPGP1 number 4786"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GFM1 number 4801"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GFPT2 number 4806"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GFRA1 number 4807"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GGA1 number 4813"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GGA2 number 4814"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GGACT number 4816"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GGCX number 4818"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GGT1 number 4822"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GGT5 number 4823"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GHITM number 4827"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GHSR number 4832"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GID4 number 4833"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GIMD1 number 4838"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GINS4 number 4844"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene GIPC3 number 4847"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GJA5 number 4855"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GJB1 number 4858"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GJB2 number 4859"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GJC1 number 4862"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GJD2 number 4864"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene GLA number 4872"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GLDN number 4878"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GLG1 number 4880"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GLI2 number 4882"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GLMN number 4888"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GLOD4 number 4890"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GLRB number 4897"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GLT8D2 number 4905"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GLYCTK number 4909"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GLYR1 number 4910"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GM2A number 4911"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GMEB2 number 4915"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GMNN number 4919"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GMPPA number 4920"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GNAT2 number 4937"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GNAT3 number 4938"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GNB2 number 4942"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GNB3 number 4943"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GNB5 number 4945"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GNG10 number 4947"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GNG3 number 4951"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene GNPTAB number 4966"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene GOLGA7 number 4975"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GOLM2 number 4980"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GOLT1B number 4984"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GOPC number 4985"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GORASP1 number 4987"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GOSR1 number 4989"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GOT1 number 4991"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GP9 number 4996"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene GPAA1 number 4998"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPATCH2L number 5008"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPATCH8 number 5011"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPBAR1 number 5012"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPM6B number 5033"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPR107 number 5039"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPR12 number 5042"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene GPR132 number 5043"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPR135 number 5044"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene GPR137 number 5045"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPR139 number 5048"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPR149 number 5054"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPR153 number 5059"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPR157 number 5062"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPR161 number 5065"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPR180 number 5074"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPR183 number 5076"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPR27 number 5083"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene GPR37L1 number 5089"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPR62 number 5097"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene GPR63 number 5098"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPR65 number 5099"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPR82 number 5103"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene GPR85 number 5106"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPR87 number 5107"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPRC5B number 5111"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPRIN3 number 5115"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPSM3 number 5119"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618392.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene GPT2 number 5120"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GPX7 number 5125"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GREM1 number 5140"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GREM2 number 5141"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GRHL1 number 5142"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GRHL2 number 5143"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GRIA1 number 5146"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GRIN1 number 5158"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GRIN3A number 5163"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GRINA number 5165"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GRPEL2 number 5185"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GRTP1 number 5188"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GRXCR2 number 5191"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GSC number 5192"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene GSDME number 5194"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GSN number 5200"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GSTO1 number 5206"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GSX1 number 5209"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene GTF2A1L number 5213"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GTF2A2 number 5214"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GTF2E1 number 5216"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GTF2H5 number 5224"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GTF3C1 number 5228"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GTF3C2 number 5229"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GTF3C3 number 5230"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GTF3C5 number 5232"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GTF3C6 number 5233"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GTPBP2 number 5236"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GTPBP6 number 5239"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GTPBP8 number 5240"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GUCA1ANB number 5242"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GUCA1B number 5243"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GUCA1C number 5244"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GXYLT1 number 5255"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene GZF1 number 5262"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HACD1 number 5268"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HACD4 number 5271"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HACL1 number 5273"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HADHA number 5275"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HAO1 number 5281"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HAO2 number 5282"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HAPLN1 number 5284"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HARBI1 number 5288"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene HAT1 number 5293"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HDGF number 5327"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HDHD5 number 5332"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HEATR3 number 5335"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HEATR4 number 5336"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HEATR5B number 5338"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HEATR6 number 5339"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HEATR9 number 5340"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HECTD3 number 5346"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HELB number 5351"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HELQ number 5353"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HELZ2 number 5356"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HENMT1 number 5359"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HEXIM1 number 5378"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene HGF number 5383"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HGSNAT number 5387"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HHLA2 number 5394"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HIBADH number 5395"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HIBCH number 5396"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HIC1 number 5397"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene HIGD1C number 5405"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene HIPK2 number 5414"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HIVEP2 number 5419"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HK1 number 5422"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HKDC1 number 5425"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HMCN1 number 5434"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HMCN2 number 5435"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HMG20A number 5436"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HMGA2 number 5439"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HMGCLL1 number 5444"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HNF4G number 5462"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HNRNPAB number 5468"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HNRNPH1 number 5472"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HNRNPLL number 5476"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HOGA1 number 5483"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HOMER1 number 5484"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HOXA11 number 5496"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene HOXA6 number 5502"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene HOXB2 number 5507"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene HOXB9 number 5514"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HOXC11 number 5516"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene HOXC4 number 5519"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene HOXC5 number 5520"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene HOXC8 number 5522"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene HOXC9 number 5523"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HPCAL1 number 5535"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HPCAL4 number 5536"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HPGD number 5539"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HPN number 5541"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HPS1 number 5543"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HPS6 number 5546"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene HPSE number 5547"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HR number 5550"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HRAS number 5551"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HS3ST2 number 5560"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HS6ST2 number 5565"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HSD17B1 number 5573"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HSDL2 number 5583"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HSPA12B number 5592"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HSPA5 number 5598"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HSPA9 number 5600"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HSPB7 number 5606"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HTR1B number 5617"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HTR1E number 5619"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene HTR6 number 5628"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HTRA1 number 5630"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HUNK number 5635"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HVCN1 number 5638"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene HYKK number 5644"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IARS1 number 5649"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene IBA57 number 5651"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ICE1 number 5656"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IDUA number 5671"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IER5L number 5676"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene IFITM5 number 5684"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene IFNAR1 number 5685"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IFNAR2 number 5686"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IFNGR1 number 5688"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IFT20 number 5697"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IFT22 number 5698"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IFT43 number 5700"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IFT80 number 5705"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IFT81 number 5706"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IGF2BP2 number 5716"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IGF2R number 5718"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IGFBP7 number 5726"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IGLON5 number 5731"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IGSF22 number 5735"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IGSF9B number 5741"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IKZF2 number 5749"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IL10RB number 5754"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IL13RA1 number 5761"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IL15RA number 5764"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IL17RD number 5773"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IL1RL1 number 5781"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IL21R number 5784"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IL22 number 5785"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene IL22RA1 number 5786"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene IL23R number 5789"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IL6 number 5800"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IL6ST number 5802"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IMPA1 number 5817"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IMPDH2 number 5821"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IMPG1 number 5822"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ING4 number 5833"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ING5 number 5834"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene INHBA number 5836"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene INO80B number 5841"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene INPP5B number 5849"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene INSIG1 number 5858"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene INSL5 number 5860"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene INSM1 number 5861"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene INTS10 number 5869"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene INTS4 number 5876"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IP6K3 number 5886"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IQCG number 5904"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IQGAP3 number 5911"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IRF2BP1 number 5922"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene IRF6 number 5927"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IRX2 number 5936"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene IRX4 number 5938"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene ISCA2 number 5941"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ISL1 number 5945"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ISLR number 5947"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene ISX number 5954"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ISY1 number 5955"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ITFG1 number 5958"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ITGA5 number 5967"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ITGA7 number 5969"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ITGA8 number 5970"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ITGB1BP1 number 5974"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ITGB1BP2 number 5975"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ITGB6 number 5980"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ITGB7 number 5981"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ITIH5 number 5985"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IVD number 6005"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IYD number 6007"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene IZUMO4 number 6010"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene JADE2 number 6012"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene JAGN1 number 6016"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene JAK2 number 6018"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene JAK3 number 6019"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene JAM2 number 6023"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene JMJD6 number 6034"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene JMY number 6037"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KANSL1 number 6056"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KAT6B number 6067"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KAT7 number 6068"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KCMF1 number 6085"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KCNA1 number 6086"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene KCNE3 number 6104"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KCNE5 number 6106"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene KCNF1 number 6107"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene KCNH3 number 6114"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KCNJ1 number 6124"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KCNJ15 number 6128"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KCNJ4 number 6132"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KCNJ5 number 6133"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KCNJ9 number 6136"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KCNK6 number 6149"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KCNN1 number 6156"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KCNU1 number 6171"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KCNV2 number 6173"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene KCTD2 number 6184"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KCTD20 number 6185"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KCTD3 number 6187"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KDM1B number 6198"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KDM3A number 6201"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KDM3B number 6202"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KDM8 number 6211"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KDR number 6212"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KIAA0232 number 6225"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KIAA0895 number 6232"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KIAA1109 number 6235"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KIAA1143 number 6236"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KIAA1522 number 6241"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KIAA1755 number 6246"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KIAA2026 number 6251"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KIF13B number 6256"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KIF19 number 6263"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KIF26A number 6274"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KIF3A number 6280"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KIFBP number 6291"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KIFC2 number 6293"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KIRREL1 number 6296"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KITLG number 6302"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KLC1 number 6306"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KLF2 number 6317"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene KLF3 number 6318"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KLF7 number 6322"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KLHDC1 number 6325"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KLHDC2 number 6327"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KLHDC3 number 6328"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KLHL13 number 6338"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KLHL14 number 6339"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KLHL26 number 6350"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KMT2E number 6376"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KNL1 number 6382"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene KPNA6 number 6391"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KREMEN1 number 6397"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KREMEN2 number 6398"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KRIT1 number 6400"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene KRT222 number 6403"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LACTB2 number 6425"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LAGE3 number 6429"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LAMA2 number 6431"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LAMB3 number 6437"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LAMB4 number 6438"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LAMC1 number 6439"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LAMC2 number 6440"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LAMC3 number 6441"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LAMP1 number 6442"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LAMP3 number 6444"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LAMP5 number 6445"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LAPTM5 number 6457"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LARP1B number 6461"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LAT2 number 6471"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LCK number 6482"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LCP2 number 6489"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LDAH number 6492"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LDB2 number 6494"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LDHB number 6497"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LDLR number 6499"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LDLRAD2 number 6501"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LEMD1 number 6509"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LENG8 number 6513"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LEPR number 6516"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LGALS8 number 6530"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LGALS9 number 6531"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LGR4 number 6537"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LGSN number 6540"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LHFPL2 number 6542"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LHX2 number 6548"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LHX3 number 6549"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LHX8 number 6553"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LIFR number 6556"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LIMA1 number 6561"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LIMK2 number 6567"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LIMS2 number 6569"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LIN7A number 6575"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LIN7C number 6577"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LIN9 number 6578"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LINGO4 number 6582"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LKAAEAR1 number 6593"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LLGL1 number 6594"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LLGL2 number 6595"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LLPH number 6596"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LMAN1 number 6597"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LMLN number 6607"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LMNB1 number 6609"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LMNB2 number 6610"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LMNTD2 number 6612"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LMO2 number 6614"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LMOD2 number 6619"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LNX1 number 6627"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LNX2 number 6628"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102929294 number 6629"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102929298 number 6631"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102929329 number 6639"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102929444 number 6661"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102929519 number 6670"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102929559 number 6676"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102929605 number 6686"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102929613 number 6687"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102929631 number 6689"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102929645 number 6693"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102929653 number 6695"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102929814 number 6717"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102929815 number 6718"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102929832 number 6724"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930039 number 6755"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930080 number 6761"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930167 number 6772"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930181 number 6778"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930214 number 6782"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC102930327 number 6803"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930346 number 6810"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930458 number 6829"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930497 number 6840"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930502 number 6842"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930516 number 6843"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930517 number 6844"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930533 number 6847"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930604 number 6856"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102930609 number 6858"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102930615 number 6860"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930625 number 6863"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930647 number 6866"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930688 number 6875"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930697 number 6877"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102930698 number 6878"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102930701 number 6879"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930765 number 6892"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930790 number 6896"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930792 number 6897"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930812 number 6900"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930830 number 6902"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102930835 number 6903"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102930909 number 6916"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102930929 number 6922"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102930967 number 6929"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931016 number 6938"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102931048 number 6947"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102931053 number 6948"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102931073 number 6955"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931081 number 6957"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931111 number 6958"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931141 number 6962"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931149 number 6964"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931185 number 6967"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102931277 number 6985"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102931294 number 6987"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931368 number 7002"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931395 number 7008"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931425 number 7015"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931536 number 7029"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102931639 number 7046"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931648 number 7047"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931701 number 7054"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102931711 number 7056"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931719 number 7059"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931746 number 7063"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931747 number 7064"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC102931790 number 7072"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931801 number 7074"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931805 number 7076"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931836 number 7082"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931859 number 7086"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931862 number 7087"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931863 number 7088"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931963 number 7105"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931965 number 7106"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102931998 number 7109"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932022 number 7111"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932027 number 7114"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932055 number 7121"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932122 number 7131"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932126 number 7132"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932174 number 7140"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932184 number 7143"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932237 number 7152"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102932238 number 7153"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102932250 number 7156"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932257 number 7159"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102932297 number 7166"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932311 number 7168"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102932337 number 7174"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932346 number 7178"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932409 number 7188"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932459 number 7194"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932523 number 7209"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102932558 number 7214"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932562 number 7216"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102932601 number 7221"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932612 number 7224"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932629 number 7226"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932677 number 7229"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932697 number 7235"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102932780 number 7249"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932805 number 7254"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932811 number 7255"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932837 number 7262"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102932929 number 7281"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102932988 number 7290"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102932998 number 7291"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102933026 number 7297"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102933055 number 7300"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102933084 number 7305"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102933104 number 7307"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102933112 number 7308"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102933117 number 7310"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102933126 number 7311"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102933224 number 7328"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102933267 number 7331"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102933309 number 7343"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102933321 number 7345"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102933374 number 7352"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102933397 number 7358"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102933496 number 7377"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102933505 number 7379"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102933579 number 7392"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102933620 number 7396"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102933627 number 7399"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102933687 number 7405"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102933707 number 7409"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102933710 number 7410"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102933733 number 7414"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102933768 number 7421"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102933796 number 7425"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102933929 number 7448"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102933975 number 7456"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102933989 number 7458"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934017 number 7467"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934045 number 7476"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934068 number 7480"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934134 number 7491"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102934150 number 7494"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934154 number 7496"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934221 number 7506"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102934223 number 7507"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934226 number 7508"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102934244 number 7512"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934268 number 7515"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934296 number 7523"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934304 number 7524"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934360 number 7536"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102934370 number 7539"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934412 number 7545"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102934596 number 7577"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934640 number 7580"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934685 number 7587"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934722 number 7593"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102934807 number 7606"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934855 number 7620"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934896 number 7628"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934908 number 7629"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934954 number 7634"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934993 number 7642"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102934999 number 7643"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935035 number 7647"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935074 number 7652"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935087 number 7654"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935115 number 7659"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102935129 number 7663"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102935130 number 7664"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC102935133 number 7665"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102935276 number 7680"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935280 number 7682"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102935328 number 7688"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935338 number 7690"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935352 number 7693"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618392.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102935353 number 7694"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935415 number 7707"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935425 number 7711"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935427 number 7712"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935443 number 7715"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935475 number 7721"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935521 number 7731"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102935545 number 7737"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935609 number 7750"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102935616 number 7752"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935635 number 7756"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102935655 number 7760"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935666 number 7764"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935673 number 7767"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935703 number 7773"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102935704 number 7774"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935738 number 7779"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102935762 number 7782"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935816 number 7796"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102935859 number 7803"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102935872 number 7809"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102935949 number 7818"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618350.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102935961 number 7819"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102935973 number 7821"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936018 number 7826"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102936042 number 7830"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936099 number 7841"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102936111 number 7844"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618379.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102936138 number 7848"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936143 number 7850"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936290 number 7871"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102936380 number 7890"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936385 number 7891"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936460 number 7903"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102936508 number 7916"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102936514 number 7917"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936516 number 7918"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102936539 number 7923"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936556 number 7926"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936557 number 7927"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936608 number 7933"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936613 number 7934"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102936625 number 7939"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102936629 number 7941"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936669 number 7947"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936721 number 7954"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936745 number 7960"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936754 number 7963"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936834 number 7979"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102936849 number 7982"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936856 number 7983"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936860 number 7985"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936895 number 7988"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102936991 number 8005"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937094 number 8025"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937124 number 8031"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102937142 number 8035"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102937143 number 8036"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937153 number 8040"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937154 number 8041"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102937166 number 8043"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937186 number 8046"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937262 number 8060"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937275 number 8062"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937323 number 8073"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937329 number 8074"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937338 number 8076"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937340 number 8078"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937342 number 8079"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937354 number 8081"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102937418 number 8089"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937440 number 8095"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102937471 number 8099"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102937485 number 8102"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618351.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102937488 number 8103"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937516 number 8108"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937518 number 8109"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937536 number 8113"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937619 number 8129"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937641 number 8135"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937661 number 8137"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937705 number 8145"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937715 number 8146"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102937734 number 8147"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937739 number 8148"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937746 number 8151"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937789 number 8162"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937840 number 8168"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937867 number 8174"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937870 number 8175"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937872 number 8176"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102937966 number 8195"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102937979 number 8198"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102938005 number 8200"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938129 number 8214"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938172 number 8222"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938230 number 8232"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102938233 number 8233"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938308 number 8244"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102938434 number 8263"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938452 number 8266"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102938512 number 8273"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938557 number 8278"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938581 number 8283"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102938586 number 8284"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938596 number 8285"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938676 number 8301"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102938683 number 8304"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938691 number 8307"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938707 number 8309"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102938714 number 8312"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938716 number 8313"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938721 number 8315"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938725 number 8316"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938736 number 8319"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938742 number 8321"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938791 number 8326"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938818 number 8331"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938875 number 8340"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938904 number 8347"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938910 number 8348"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938913 number 8349"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC102938920 number 8352"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938936 number 8357"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938945 number 8361"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102938962 number 8368"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102939028 number 8379"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102939044 number 8383"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939094 number 8387"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939112 number 8390"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102939129 number 8395"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939179 number 8405"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939195 number 8410"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939243 number 8418"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102939257 number 8422"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939266 number 8424"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939271 number 8426"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939310 number 8429"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939324 number 8433"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102939330 number 8434"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939338 number 8436"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939346 number 8439"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102939352 number 8440"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939462 number 8459"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102939469 number 8461"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC102939482 number 8466"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102939494 number 8468"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939502 number 8469"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC102939531 number 8473"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939620 number 8482"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939664 number 8490"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939665 number 8491"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939692 number 8500"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939696 number 8501"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102939778 number 8510"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102939782 number 8512"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939851 number 8526"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102939869 number 8527"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939895 number 8530"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102939922 number 8537"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939940 number 8541"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102939959 number 8545"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102939967 number 8547"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939969 number 8548"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939970 number 8549"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939983 number 8552"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102939992 number 8554"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940046 number 8562"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102940092 number 8572"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940134 number 8580"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102940161 number 8586"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102940179 number 8589"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940187 number 8590"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940192 number 8592"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940202 number 8595"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940268 number 8607"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940278 number 8611"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940281 number 8612"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940329 number 8618"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940341 number 8620"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940420 number 8632"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102940422 number 8633"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940455 number 8636"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940472 number 8639"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940473 number 8640"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940475 number 8641"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940499 number 8648"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940517 number 8651"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102940526 number 8653"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940584 number 8666"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940592 number 8667"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940605 number 8668"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940606 number 8669"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102940647 number 8675"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940649 number 8676"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102940652 number 8677"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940657 number 8679"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102940697 number 8683"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940714 number 8689"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940718 number 8690"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940738 number 8693"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940789 number 8706"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940797 number 8708"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102940858 number 8712"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940914 number 8721"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940929 number 8726"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102940981 number 8734"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102940984 number 8735"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941010 number 8739"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941079 number 8753"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941199 number 8767"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941203 number 8769"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941213 number 8771"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941219 number 8773"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941257 number 8780"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102941277 number 8783"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941295 number 8787"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102941337 number 8793"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941343 number 8796"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941371 number 8804"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC102941380 number 8807"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941397 number 8810"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102941417 number 8813"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941435 number 8818"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102941452 number 8824"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102941511 number 8836"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941514 number 8838"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102941516 number 8839"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102941531 number 8842"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102941574 number 8848"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941576 number 8849"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941598 number 8851"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941612 number 8852"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941632 number 8854"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941655 number 8859"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941682 number 8863"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941700 number 8868"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941765 number 8884"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618385.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102941785 number 8888"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102941811 number 8890"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941816 number 8893"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941826 number 8895"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102941838 number 8898"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941857 number 8901"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941873 number 8905"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941881 number 8907"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941909 number 8910"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941939 number 8914"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102941946 number 8915"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941965 number 8920"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941968 number 8921"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102941991 number 8925"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102942038 number 8931"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942059 number 8936"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942063 number 8938"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942127 number 8950"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942136 number 8952"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942139 number 8953"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102942208 number 8964"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942209 number 8965"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942285 number 8975"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942325 number 8977"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942326 number 8978"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942332 number 8980"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102942355 number 8991"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618375.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102942396 number 8996"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942398 number 8997"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942400 number 8998"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942405 number 8999"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942503 number 9014"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102942532 number 9016"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942534 number 9017"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942577 number 9026"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942580 number 9027"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942655 number 9034"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942671 number 9037"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102942672 number 9038"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102942673 number 9039"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942686 number 9041"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942771 number 9056"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102942773 number 9057"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942794 number 9061"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942795 number 9062"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942801 number 9063"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102942858 number 9073"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942881 number 9080"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102942939 number 9087"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102942970 number 9090"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943028 number 9100"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943051 number 9105"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943084 number 9109"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102943120 number 9113"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943152 number 9117"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943176 number 9122"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943243 number 9135"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943252 number 9138"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943313 number 9148"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102943321 number 9151"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943322 number 9152"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC102943365 number 9159"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943367 number 9160"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943373 number 9161"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102943407 number 9167"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943445 number 9173"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943456 number 9176"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102943470 number 9181"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943586 number 9197"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943671 number 9209"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC102943686 number 9211"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC102943727 number 9217"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943741 number 9220"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943751 number 9222"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943755 number 9223"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943760 number 9225"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943777 number 9226"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943815 number 9232"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102943825 number 9237"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943847 number 9244"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943903 number 9254"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102943912 number 9255"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102943966 number 9263"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944006 number 9270"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944065 number 9276"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944092 number 9279"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944114 number 9284"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944152 number 9290"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944173 number 9292"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944268 number 9300"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102944271 number 9301"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944277 number 9302"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944278 number 9303"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102944284 number 9305"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944330 number 9309"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944361 number 9318"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944368 number 9319"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102944375 number 9321"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944426 number 9332"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102944537 number 9358"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944584 number 9366"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944634 number 9372"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102944636 number 9373"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944748 number 9390"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944784 number 9398"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944830 number 9408"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944862 number 9412"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944882 number 9419"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102944926 number 9424"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102944948 number 9433"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102944998 number 9440"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945012 number 9444"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945028 number 9452"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102945042 number 9455"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945058 number 9459"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945076 number 9461"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945093 number 9464"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945101 number 9465"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102945131 number 9469"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945151 number 9475"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102945178 number 9481"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102945196 number 9483"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102945219 number 9485"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945232 number 9487"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945298 number 9499"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945363 number 9507"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945366 number 9508"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102945393 number 9516"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102945398 number 9517"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945437 number 9522"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945518 number 9541"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102945523 number 9542"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945548 number 9546"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945576 number 9550"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945642 number 9564"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102945644 number 9566"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945654 number 9567"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102945655 number 9568"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945694 number 9577"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945695 number 9578"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945757 number 9588"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945807 number 9595"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102945845 number 9605"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945900 number 9615"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945926 number 9620"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945929 number 9621"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945930 number 9622"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102945948 number 9625"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102945962 number 9628"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102945967 number 9630"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102945995 number 9636"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102946053 number 9642"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946063 number 9646"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618370.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102946105 number 9653"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946106 number 9654"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946167 number 9672"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946201 number 9677"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946237 number 9683"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102946280 number 9692"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102946281 number 9693"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946308 number 9702"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946316 number 9703"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946326 number 9705"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946339 number 9710"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946353 number 9711"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102946478 number 9727"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946489 number 9729"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946497 number 9733"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102946526 number 9743"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102946545 number 9747"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC102946588 number 9752"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102946593 number 9755"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946613 number 9760"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102946615 number 9761"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946640 number 9763"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102946652 number 9766"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102946657 number 9770"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102946737 number 9783"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102946742 number 9787"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946787 number 9792"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946818 number 9798"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946884 number 9807"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946904 number 9813"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946939 number 9818"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946959 number 9822"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946960 number 9823"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102946980 number 9829"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102946998 number 9833"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947012 number 9837"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947047 number 9840"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102947048 number 9841"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102947095 number 9850"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947119 number 9854"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947120 number 9855"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947137 number 9858"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947141 number 9859"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947196 number 9869"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947268 number 9879"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102947274 number 9882"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947282 number 9883"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947314 number 9892"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102947325 number 9893"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947327 number 9894"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947375 number 9903"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102947484 number 9924"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947489 number 9926"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947510 number 9928"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947518 number 9931"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947623 number 9948"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102947650 number 9953"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947665 number 9955"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947668 number 9956"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947706 number 9965"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947707 number 9966"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947745 number 9977"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947761 number 9979"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102947763 number 9980"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947766 number 9981"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947869 number 10007"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC102947883 number 10010"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102947951 number 10026"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102947970 number 10028"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102947981 number 10029"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102948026 number 10037"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102948041 number 10040"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102948081 number 10046"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC102948100 number 10050"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102948104 number 10052"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102948129 number 10056"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102948136 number 10058"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102948146 number 10060"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102948152 number 10061"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618375.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102948174 number 10068"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102948216 number 10077"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102948218 number 10078"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102948226 number 10079"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102948253 number 10086"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC102948298 number 10093"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102948335 number 10103"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC102948381 number 10111"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018199 number 10117"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018202 number 10118"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018211 number 10122"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018212 number 10123"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018215 number 10125"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018225 number 10132"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018261 number 10144"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018275 number 10153"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018280 number 10158"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018314 number 10176"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018336 number 10178"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018342 number 10181"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018345 number 10182"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018351 number 10185"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018353 number 10187"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC114018362 number 10191"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018384 number 10204"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018391 number 10205"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018393 number 10206"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018401 number 10209"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018449 number 10226"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018450 number 10227"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018475 number 10231"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018486 number 10234"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018499 number 10241"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018508 number 10244"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018512 number 10247"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018517 number 10251"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018530 number 10256"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018540 number 10261"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018546 number 10264"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018587 number 10280"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018602 number 10290"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018625 number 10296"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018646 number 10302"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018657 number 10304"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018671 number 10311"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018677 number 10314"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018695 number 10320"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018724 number 10328"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018729 number 10331"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618375.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018730 number 10332"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618375.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018750 number 10347"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018765 number 10353"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018782 number 10361"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018784 number 10362"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018787 number 10364"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018789 number 10365"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018817 number 10370"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018838 number 10374"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018866 number 10385"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018917 number 10403"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018928 number 10406"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114018935 number 10410"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018956 number 10421"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018973 number 10427"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018980 number 10429"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114018991 number 10430"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019056 number 10460"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019078 number 10465"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019088 number 10469"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019154 number 10491"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019155 number 10492"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019161 number 10494"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019173 number 10498"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019215 number 10513"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019228 number 10522"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019264 number 10539"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019268 number 10541"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019281 number 10545"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019291 number 10548"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019319 number 10558"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019327 number 10564"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019356 number 10578"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019365 number 10581"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019376 number 10585"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019380 number 10589"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019385 number 10592"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019422 number 10605"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019456 number 10623"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019460 number 10626"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019502 number 10642"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019513 number 10646"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019537 number 10656"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019551 number 10662"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019552 number 10663"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019570 number 10669"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019587 number 10674"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019598 number 10678"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019602 number 10680"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019605 number 10681"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019654 number 10697"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019703 number 10706"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019717 number 10711"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019726 number 10715"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019738 number 10722"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019764 number 10737"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019774 number 10743"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019780 number 10744"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019789 number 10745"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019811 number 10757"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019849 number 10770"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019855 number 10773"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019870 number 10774"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019871 number 10775"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019908 number 10793"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019919 number 10794"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618392.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019926 number 10796"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019954 number 10811"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019963 number 10815"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114019974 number 10817"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019978 number 10821"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114019990 number 10826"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114020020 number 10839"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020041 number 10842"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114020050 number 10849"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020073 number 10858"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020076 number 10860"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020114 number 10877"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020116 number 10879"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020138 number 10890"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020169 number 10906"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020206 number 10922"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114020251 number 10943"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020254 number 10946"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020343 number 10972"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020363 number 10978"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114020379 number 10986"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114020393 number 10989"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020397 number 10991"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020402 number 10993"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020404 number 10994"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020408 number 10997"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020409 number 10998"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020410 number 10999"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020435 number 11010"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020495 number 11034"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020516 number 11041"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020518 number 11042"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020522 number 11046"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020558 number 11061"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114020565 number 11064"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114020584 number 11077"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020621 number 11087"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020651 number 11101"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020671 number 11108"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020672 number 11109"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020717 number 11116"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020719 number 11117"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020724 number 11119"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020748 number 11127"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020749 number 11128"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020758 number 11131"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020798 number 11147"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114020815 number 11154"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020829 number 11159"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114020851 number 11167"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020852 number 11168"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020879 number 11173"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114020891 number 11178"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114020894 number 11179"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020912 number 11187"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114020919 number 11191"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114020925 number 11195"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020966 number 11213"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114020978 number 11215"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021006 number 11225"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021007 number 11226"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021015 number 11231"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021047 number 11242"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021061 number 11247"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021082 number 11256"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021086 number 11258"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021100 number 11267"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021103 number 11268"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021154 number 11292"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021168 number 11298"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021186 number 11308"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021219 number 11316"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021238 number 11320"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021270 number 11331"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021276 number 11334"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021277 number 11335"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021282 number 11337"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021291 number 11341"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021304 number 11346"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021305 number 11347"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021315 number 11352"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021353 number 11367"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021383 number 11380"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021388 number 11381"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021408 number 11390"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021418 number 11395"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021433 number 11400"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021444 number 11401"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021460 number 11405"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021466 number 11408"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021505 number 11422"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021517 number 11426"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021522 number 11427"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021532 number 11432"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021559 number 11441"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021570 number 11445"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021580 number 11449"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021596 number 11455"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021600 number 11457"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021602 number 11459"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021652 number 11477"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021670 number 11485"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021686 number 11491"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021689 number 11494"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021714 number 11512"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021715 number 11513"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021723 number 11515"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021759 number 11522"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021846 number 11553"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021849 number 11554"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021852 number 11555"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021885 number 11561"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021893 number 11565"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021898 number 11568"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021901 number 11570"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021918 number 11583"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021930 number 11588"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021932 number 11589"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021946 number 11591"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021958 number 11596"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114021969 number 11604"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021972 number 11606"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021986 number 11608"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114021996 number 11613"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022024 number 11630"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022028 number 11632"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114022041 number 11634"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114022044 number 11635"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022047 number 11636"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022052 number 11639"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022058 number 11642"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022064 number 11646"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022086 number 11648"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022102 number 11655"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114022116 number 11661"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114022129 number 11665"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114022142 number 11671"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022152 number 11675"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022166 number 11680"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114022198 number 11688"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022199 number 11689"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022223 number 11699"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022239 number 11706"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114022241 number 11707"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022245 number 11709"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022257 number 11716"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022264 number 11718"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114022268 number 11720"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022295 number 11725"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022300 number 11728"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022304 number 11730"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC114022387 number 11766"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114022390 number 11767"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022399 number 11771"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022414 number 11782"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022415 number 11783"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022425 number 11788"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022431 number 11790"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022444 number 11793"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC114022468 number 11802"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022492 number 11807"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022511 number 11817"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022543 number 11829"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022545 number 11831"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022565 number 11845"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022589 number 11857"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022590 number 11858"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC114022601 number 11860"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022609 number 11864"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022610 number 11865"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022627 number 11869"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022636 number 11874"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022664 number 11883"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022667 number 11884"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022679 number 11887"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC114022684 number 11889"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563548 number 11896"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563552 number 11899"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563553 number 11900"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563571 number 11912"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563575 number 11916"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563578 number 11918"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563581 number 11921"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563588 number 11926"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563592 number 11929"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563594 number 11931"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119563604 number 11941"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119563605 number 11942"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563608 number 11945"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563609 number 11946"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563615 number 11952"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563623 number 11960"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563628 number 11963"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563629 number 11964"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563648 number 11971"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563652 number 11975"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563668 number 11985"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563681 number 11997"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563689 number 12005"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563692 number 12008"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563701 number 12014"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563703 number 12016"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563720 number 12029"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563723 number 12032"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563737 number 12042"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563758 number 12046"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563774 number 12062"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563787 number 12074"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563816 number 12100"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563818 number 12102"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563825 number 12109"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563837 number 12118"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563838 number 12119"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563844 number 12125"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563849 number 12130"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563851 number 12132"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563855 number 12136"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563858 number 12138"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563874 number 12154"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563881 number 12161"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563883 number 12163"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563888 number 12168"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563891 number 12171"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563893 number 12172"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563897 number 12176"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563898 number 12177"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563899 number 12178"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563901 number 12180"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563907 number 12184"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563911 number 12187"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563921 number 12193"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563924 number 12196"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563929 number 12199"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563934 number 12204"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563937 number 12207"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563945 number 12214"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563947 number 12216"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563955 number 12223"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563960 number 12227"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563965 number 12231"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563969 number 12235"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563976 number 12240"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119563979 number 12243"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563983 number 12247"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119563988 number 12252"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564006 number 12267"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564015 number 12276"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564061 number 12305"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564062 number 12306"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564063 number 12307"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564067 number 12311"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564072 number 12316"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564077 number 12321"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564079 number 12323"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564082 number 12326"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564085 number 12329"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564095 number 12338"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564112 number 12354"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564117 number 12358"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564118 number 12359"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564133 number 12374"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564136 number 12377"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119564156 number 12395"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564159 number 12397"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564163 number 12401"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564165 number 12402"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564169 number 12406"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564173 number 12409"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564174 number 12410"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564176 number 12412"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564184 number 12418"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564186 number 12419"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564192 number 12423"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564210 number 12438"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564215 number 12443"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564216 number 12444"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564218 number 12445"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564219 number 12446"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564222 number 12449"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564238 number 12462"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564241 number 12465"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119564242 number 12466"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119564244 number 12467"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564257 number 12479"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564281 number 12503"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564293 number 12512"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564296 number 12514"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564307 number 12522"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564316 number 12527"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564317 number 12528"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564319 number 12530"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564323 number 12532"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564330 number 12538"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564342 number 12548"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564352 number 12554"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564357 number 12557"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564360 number 12559"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564361 number 12560"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564387 number 12576"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564401 number 12588"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564407 number 12594"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564420 number 12607"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564424 number 12608"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564427 number 12611"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564431 number 12615"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564441 number 12624"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564447 number 12629"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564455 number 12635"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564458 number 12638"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564463 number 12643"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564465 number 12644"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564466 number 12645"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564469 number 12647"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564476 number 12652"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564480 number 12656"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564489 number 12665"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564503 number 12676"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564507 number 12680"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564517 number 12689"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564519 number 12691"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564523 number 12695"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564525 number 12697"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564529 number 12700"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564533 number 12702"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564534 number 12703"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564541 number 12710"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564542 number 12711"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564544 number 12713"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564555 number 12721"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564562 number 12728"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564565 number 12730"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564602 number 12758"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564606 number 12762"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119564609 number 12763"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564619 number 12773"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564623 number 12777"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564632 number 12782"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564636 number 12786"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564658 number 12799"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564721 number 12803"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564757 number 12804"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564787 number 12808"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564797 number 12816"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564827 number 12842"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564841 number 12852"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564844 number 12855"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564866 number 12868"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564874 number 12874"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564881 number 12881"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119564884 number 12883"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564894 number 12893"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564896 number 12894"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564910 number 12902"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564917 number 12908"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119564918 number 12909"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564925 number 12914"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564933 number 12916"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564952 number 12930"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564960 number 12934"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618342.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564974 number 12942"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618357.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564975 number 12943"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618358.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564977 number 12945"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564978 number 12946"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564983 number 12951"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618370.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564986 number 12953"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564989 number 12955"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119564990 number 12956"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119564993 number 12959"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565000 number 12965"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565003 number 12966"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565006 number 12967"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565012 number 12969"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618375.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565015 number 12971"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565016 number 12972"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565029 number 12976"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565045 number 12986"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565047 number 12988"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618391.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565054 number 12990"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565056 number 12991"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618400.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565057 number 12992"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618401.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565058 number 12993"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565072 number 13005"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565076 number 13009"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565080 number 13013"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565083 number 13016"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565087 number 13020"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565122 number 13049"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565128 number 13055"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565131 number 13058"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565138 number 13064"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119565159 number 13084"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565160 number 13085"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565161 number 13086"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119565168 number 13092"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565197 number 13115"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565206 number 13124"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565210 number 13128"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565211 number 13129"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565214 number 13132"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565215 number 13133"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565216 number 13134"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565222 number 13139"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565228 number 13145"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565229 number 13146"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565232 number 13149"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565240 number 13157"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565242 number 13159"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565243 number 13160"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565253 number 13170"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565258 number 13174"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565267 number 13180"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565291 number 13203"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565297 number 13209"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565310 number 13221"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565322 number 13231"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565350 number 13239"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565355 number 13244"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565356 number 13245"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565357 number 13246"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565358 number 13247"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565360 number 13249"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565362 number 13251"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565367 number 13256"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565368 number 13257"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565369 number 13258"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565374 number 13263"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565380 number 13267"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565391 number 13278"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565394 number 13281"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565398 number 13285"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119565400 number 13287"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565402 number 13289"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565407 number 13293"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565423 number 13309"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565431 number 13315"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565439 number 13322"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565446 number 13329"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565454 number 13337"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565468 number 13346"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565469 number 13347"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565472 number 13350"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565479 number 13354"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565482 number 13357"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565488 number 13360"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565496 number 13366"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565507 number 13377"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565524 number 13392"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565528 number 13396"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565529 number 13397"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565533 number 13401"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119565535 number 13403"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565542 number 13409"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565546 number 13413"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565547 number 13414"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565578 number 13440"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565584 number 13445"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565588 number 13448"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565592 number 13452"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565599 number 13458"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565601 number 13459"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565605 number 13463"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565618 number 13474"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565643 number 13496"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565677 number 13517"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565682 number 13522"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565683 number 13523"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565687 number 13527"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565696 number 13535"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565702 number 13541"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565715 number 13551"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565718 number 13554"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565719 number 13555"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565731 number 13565"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565751 number 13583"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565752 number 13584"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565761 number 13593"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565762 number 13594"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565769 number 13601"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565781 number 13613"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565784 number 13616"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565795 number 13626"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565799 number 13629"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565811 number 13640"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565814 number 13642"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565816 number 13644"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565817 number 13645"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565834 number 13658"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565835 number 13659"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565847 number 13669"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565850 number 13672"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565853 number 13675"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565864 number 13686"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565976 number 13709"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565978 number 13711"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119565988 number 13721"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565993 number 13726"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119565999 number 13732"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566041 number 13769"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566052 number 13780"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566053 number 13781"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566054 number 13782"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566056 number 13784"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566061 number 13788"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566064 number 13791"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566067 number 13794"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566068 number 13795"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566088 number 13813"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566096 number 13821"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566103 number 13828"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566108 number 13832"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566113 number 13836"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566115 number 13838"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566130 number 13851"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566134 number 13854"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566142 number 13862"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566145 number 13865"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566146 number 13866"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566164 number 13883"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566171 number 13888"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566175 number 13892"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566179 number 13896"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566208 number 13900"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566215 number 13906"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566217 number 13908"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566220 number 13911"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566231 number 13922"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566234 number 13925"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119566252 number 13939"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566256 number 13942"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119566272 number 13952"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566279 number 13956"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119566303 number 13968"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566304 number 13969"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566305 number 13970"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566309 number 13972"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566318 number 13978"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566322 number 13981"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566327 number 13986"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566328 number 13987"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566329 number 13988"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566330 number 13989"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119566353 number 14005"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566356 number 14008"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566370 number 14018"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119566375 number 14023"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566377 number 14025"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566380 number 14028"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566381 number 14029"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566383 number 14031"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566396 number 14044"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566398 number 14046"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566400 number 14048"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566408 number 14056"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566409 number 14057"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566414 number 14061"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566428 number 14073"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566430 number 14075"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566453 number 14094"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566457 number 14098"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566463 number 14103"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566473 number 14112"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566474 number 14113"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566481 number 14119"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566488 number 14124"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566490 number 14126"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566491 number 14127"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566492 number 14128"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566498 number 14131"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566503 number 14134"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566504 number 14135"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119566506 number 14137"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566518 number 14146"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566524 number 14150"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566535 number 14161"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566546 number 14171"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566548 number 14173"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566550 number 14175"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566556 number 14181"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566569 number 14193"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566575 number 14199"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566577 number 14201"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566578 number 14202"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566581 number 14204"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566590 number 14213"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566601 number 14223"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119566604 number 14225"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566613 number 14233"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566614 number 14234"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566632 number 14248"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566633 number 14249"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566652 number 14267"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566654 number 14269"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566656 number 14271"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566663 number 14277"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566665 number 14279"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566675 number 14288"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566676 number 14289"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566712 number 14325"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566713 number 14326"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119566715 number 14328"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566737 number 14345"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566740 number 14348"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566754 number 14360"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566755 number 14361"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566760 number 14366"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566763 number 14369"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566765 number 14371"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566773 number 14378"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566787 number 14390"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566792 number 14395"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566794 number 14397"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119566795 number 14398"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566799 number 14402"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566800 number 14403"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566801 number 14404"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566805 number 14408"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566810 number 14413"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566812 number 14415"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566818 number 14421"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566822 number 14424"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566824 number 14426"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566833 number 14434"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566843 number 14442"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566848 number 14446"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566852 number 14450"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566858 number 14455"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566869 number 14465"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566885 number 14480"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566888 number 14483"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566896 number 14489"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566908 number 14501"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566910 number 14503"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566911 number 14504"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566912 number 14505"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566913 number 14506"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566915 number 14508"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566941 number 14529"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566955 number 14540"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566961 number 14546"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566967 number 14552"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119566969 number 14554"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566973 number 14557"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566977 number 14559"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119566989 number 14570"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567007 number 14584"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567015 number 14591"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567019 number 14595"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567034 number 14608"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567036 number 14610"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567053 number 14627"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567064 number 14638"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567069 number 14643"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567072 number 14646"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567077 number 14650"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567079 number 14652"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567091 number 14662"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567100 number 14669"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567109 number 14678"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567114 number 14683"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567117 number 14686"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567119 number 14688"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567126 number 14692"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567127 number 14693"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567132 number 14697"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567136 number 14701"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567138 number 14703"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567140 number 14705"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567141 number 14706"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567145 number 14709"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567153 number 14716"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567171 number 14731"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567172 number 14732"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567176 number 14736"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567186 number 14745"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567187 number 14746"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567196 number 14755"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567200 number 14758"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567203 number 14761"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567204 number 14762"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567218 number 14774"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567233 number 14788"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567234 number 14789"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567242 number 14795"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567267 number 14818"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567277 number 14827"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567280 number 14830"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567284 number 14833"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567309 number 14854"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567310 number 14855"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567315 number 14860"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567323 number 14868"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567326 number 14871"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567334 number 14879"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567335 number 14880"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567340 number 14885"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567341 number 14886"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567348 number 14892"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567349 number 14893"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567355 number 14897"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119567362 number 14904"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567364 number 14906"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567365 number 14907"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567367 number 14909"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567368 number 14910"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567371 number 14913"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567378 number 14920"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567379 number 14921"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567382 number 14924"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119567384 number 14925"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567394 number 14933"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567399 number 14937"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567402 number 14939"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567403 number 14940"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567404 number 14941"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567405 number 14942"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567414 number 14950"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567417 number 14953"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567427 number 14963"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567429 number 14965"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567453 number 14980"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567454 number 14981"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567456 number 14983"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567465 number 14992"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567468 number 14995"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567480 number 15005"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119567484 number 15009"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567485 number 15010"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567490 number 15015"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567507 number 15031"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567508 number 15032"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567510 number 15034"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567511 number 15035"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567518 number 15042"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567523 number 15047"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567543 number 15060"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567550 number 15067"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567560 number 15076"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567566 number 15082"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567570 number 15085"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567573 number 15088"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567592 number 15102"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567608 number 15118"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567609 number 15119"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567610 number 15120"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567616 number 15126"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567622 number 15132"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567624 number 15134"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567627 number 15137"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567643 number 15150"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567665 number 15167"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567668 number 15168"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567675 number 15173"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119567676 number 15174"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567692 number 15180"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567697 number 15185"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119567709 number 15196"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567713 number 15199"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567721 number 15205"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567726 number 15208"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567730 number 15212"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567735 number 15217"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567744 number 15226"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567745 number 15227"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567746 number 15228"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567752 number 15234"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567754 number 15236"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567757 number 15239"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567759 number 15241"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567762 number 15244"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567766 number 15248"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567779 number 15261"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567780 number 15262"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567783 number 15265"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567787 number 15269"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567791 number 15273"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567799 number 15276"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567810 number 15287"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567821 number 15294"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567822 number 15295"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567829 number 15301"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567836 number 15307"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567837 number 15308"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567842 number 15313"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567851 number 15319"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567859 number 15324"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567872 number 15333"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567873 number 15334"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567878 number 15339"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567884 number 15345"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567891 number 15352"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567894 number 15355"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567900 number 15361"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567905 number 15366"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567906 number 15367"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567914 number 15375"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567918 number 15379"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567926 number 15385"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567929 number 15387"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567937 number 15395"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567944 number 15400"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567946 number 15402"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567952 number 15408"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567962 number 15417"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567977 number 15430"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567978 number 15431"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567980 number 15432"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567989 number 15440"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119567993 number 15443"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119567996 number 15446"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene LOC119567997 number 15447"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119567998 number 15448"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119568003 number 15453"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119568004 number 15454"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119568006 number 15456"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LOC119568007 number 15457"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LOC119568022 number 15472"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LONRF1 number 15475"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LONRF3 number 15477"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LRFN1 number 15507"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LRGUK number 15512"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LRIG3 number 15516"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LRIT2 number 15518"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LRP1 number 15521"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LRP11 number 15523"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LRP3 number 15527"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LRRC23 number 15546"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LRRC25 number 15548"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LRRC32 number 15556"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LRRC3B number 15562"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LRRC4 number 15564"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene LRRC43 number 15568"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LRRC6 number 15582"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LRRC69 number 15585"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LRRC74B number 15591"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LRRC8C number 15596"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LRRC9 number 15598"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LRRN2 number 15609"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LSM14B number 15627"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LSP1 number 15637"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LSR number 15638"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LTBP2 number 15644"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LUC7L3 number 15652"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LURAP1L number 15655"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LVRN number 15658"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LY96 number 15662"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LYPD1 number 15666"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LYRM9 number 15676"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LZTFL1 number 15685"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene LZTS1 number 15687"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAEL number 15707"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAFA number 15710"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene MAK16 number 15724"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAMDC2 number 15731"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAN1A1 number 15738"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAN2B2 number 15745"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MANSC1 number 15752"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAOA number 15754"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAP1LC3A number 15758"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAP2K3 number 15765"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAP2K6 number 15768"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAP2K7 number 15769"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAP3K13 number 15774"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAP3K14 number 15775"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAP3K21 number 15780"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAP3K3 number 15781"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAP3K5 number 15783"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAP7D1 number 15798"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAPK1 number 15802"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAPK15 number 15808"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAPRE1 number 15824"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MARCHF2 number 15831"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MARCHF6 number 15835"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MARCHF8 number 15837"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MARCKSL1 number 15840"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MARF1 number 15842"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MARK2 number 15844"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MARK3 number 15845"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MARS1 number 15847"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MARVELD2 number 15850"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MASP1 number 15852"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAST3 number 15856"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAT2A number 15860"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MAZ number 15870"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MBD5 number 15877"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MBNL1 number 15882"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MBTPS1 number 15889"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MBTPS2 number 15890"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MC1R number 15891"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene MCOLN1 number 15919"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MCOLN2 number 15920"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MDC1 number 15932"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MDH1 number 15938"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MECOM number 15953"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MED1 number 15955"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MED10 number 15956"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene MED14 number 15961"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MED4 number 15978"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MED6 number 15979"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MEDAG number 15983"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MEF2C number 15985"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MEGF10 number 15987"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MEIG1 number 15994"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MEOX1 number 16004"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene METAP1D number 16013"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene METAP2 number 16014"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene METRNL number 16016"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene METTL23 number 16027"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene METTL26 number 16030"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene MEX3A number 16038"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MEX3D number 16041"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MFAP3 number 16044"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MFSD11 number 16055"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MFSD14A number 16058"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MGA number 16070"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MGAT3 number 16073"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MGRN1 number 16083"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MIB2 number 16091"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MICAL1 number 16092"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MICOS10 number 16097"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MID1IP1 number 16103"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MIDN number 16105"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MIS12 number 16128"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MIS18A number 16129"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MIS18BP1 number 16130"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MKI67 number 16136"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MKNK1 number 16139"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MLC1 number 16147"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MLLT11 number 16155"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene MLYCD number 16164"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MMACHC number 16167"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MMP15 number 16177"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MMP17 number 16179"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MMP2 number 16181"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MOB1A number 16199"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MORN1 number 16219"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MORN3 number 16221"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MORN4 number 16222"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MOSPD1 number 16226"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MPC1 number 16231"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MPC2 number 16232"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MPDZ number 16233"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene MPHOSPH8 number 16238"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MPND number 16242"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MPP5 number 16247"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MPZL1 number 16258"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MREG number 16265"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MRPL15 number 16279"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MRPL17 number 16281"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene MRPL27 number 16290"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MRPL40 number 16303"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MRPL48 number 16311"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MRPS10 number 16322"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MRPS17 number 16328"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MRPS18C number 16331"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MRPS2 number 16332"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MRPS27 number 16338"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MRPS31 number 16340"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MRPS9 number 16348"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MRTFB number 16352"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MSANTD1 number 16355"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MSC number 16359"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene MSI1 number 16366"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MSL1 number 16368"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MSMO1 number 16372"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MSRB3 number 16378"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MST1 number 16380"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MTCP1 number 16393"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene MTERF2 number 16396"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MTERF3 number 16397"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MTFR1L number 16404"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MTMR12 number 16422"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MTMR14 number 16423"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MTR number 16436"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MTRF1L number 16440"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MTUS1 number 16445"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MUC1 number 16450"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MUC13 number 16451"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MUC3A number 16454"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MUC6 number 16457"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MXD4 number 16467"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MXI1 number 16468"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene MXRA5 number 16469"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MYBPH number 16480"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MYCT1 number 16485"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MYEF2 number 16488"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MYL12B number 16499"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MYL4 number 16502"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MYL7 number 16504"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MYLK3 number 16508"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MYLK4 number 16509"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MYO1A number 16520"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MYO3B number 16529"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MYO5C number 16532"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MYO6 number 16533"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MYO7B number 16535"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MYOG number 16541"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MYOM1 number 16542"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MYOM3 number 16544"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MZB1 number 16555"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene MZT1 number 16556"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene N4BP1 number 16557"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene N4BP2 number 16558"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NAA20 number 16565"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NAAA number 16574"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NAB2 number 16578"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NACA number 16581"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NANOS1 number 16597"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene NAPRT number 16606"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NARS1 number 16608"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NAT14 number 16611"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NAT8 number 16612"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NAV1 number 16616"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NCAPH2 number 16636"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NCKAP5L number 16651"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NCLN number 16654"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NCOA6 number 16660"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NDE1 number 16670"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NDNF number 16674"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NDRG4 number 16679"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NDUFB4 number 16711"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NDUFB8 number 16715"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NDUFS5 number 16723"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NDUFS7 number 16725"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NDUFV2 number 16728"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NDUFV3 number 16729"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NEDD4L number 16743"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NEIL2 number 16750"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NEK11 number 16754"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NELFCD number 16765"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NEMP2 number 16770"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NET1 number 16775"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NETO1 number 16776"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NEURL3 number 16783"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NEUROD6 number 16787"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NFE2L2 number 16804"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NFE2L3 number 16805"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NFKBIB number 16815"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NFKBID number 16816"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NFRKB number 16820"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NFS1 number 16821"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NFU1 number 16822"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NGEF number 16829"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NGF number 16830"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NGLY1 number 16832"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NHLRC3 number 16839"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NHP2 number 16841"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NINJ1 number 16855"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NKPD1 number 16879"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NKX2-2 number 16885"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NKX2-8 number 16890"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene NLE1 number 16896"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NLGN2 number 16898"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618375.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene NLGN3 number 16899"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NLN number 16902"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NLRC3 number 16903"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NLRX1 number 16905"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NMD3 number 16908"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NME4 number 16910"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NMI number 16916"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NMNAT2 number 16918"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NMT1 number 16922"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NMU number 16924"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NMUR2 number 16926"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NOA1 number 16928"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NOC2L number 16930"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NOL12 number 16939"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene NOL3 number 16940"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NOLC1 number 16947"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NOP16 number 16952"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NOP2 number 16953"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NOP58 number 16956"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NOS2 number 16959"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NOS3 number 16960"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NOTCH1 number 16962"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NOTCH4 number 16965"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618392.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene NOXO1 number 16974"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NPAS4 number 16979"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NPB number 16980"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene NPBWR1 number 16981"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene NPFF number 16989"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NPHP1 number 16992"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NPL number 16997"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NPLOC4 number 16998"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NPR2 number 17005"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NPS number 17009"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NPY2R number 17019"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene NR1H4 number 17028"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NR2C2 number 17031"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NRBP1 number 17047"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NRBP2 number 17048"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NRF1 number 17052"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NRM number 17062"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene NRN1 number 17063"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NRN1L number 17064"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NRXN2 number 17072"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NSL1 number 17082"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NSMCE1 number 17084"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NSMF number 17087"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NSUN2 number 17089"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NT5DC3 number 17101"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NTN1 number 17109"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NTN3 number 17110"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NTN5 number 17112"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NUAK2 number 17122"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NUCB1 number 17126"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NUDCD2 number 17131"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NUDCD3 number 17132"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NUDT12 number 17134"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NUDT14 number 17136"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NUDT18 number 17140"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NUDT2 number 17142"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NUMBL number 17158"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NUP107 number 17159"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NUP153 number 17161"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NUP210L number 17167"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NUP50 number 17173"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NUP58 number 17175"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NUP98 number 17180"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NYAP1 number 17200"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NYAP2 number 17201"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene NYX number 17202"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene OAZ2 number 17207"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OCRL number 17218"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ODAD1 number 17219"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OGDHL number 17236"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OGFOD3 number 17239"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OIT3 number 17246"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OLFML1 number 17253"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OLFML2A number 17254"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OPN3 number 17272"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OPRD1 number 17275"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OPRL1 number 17277"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OPTC number 17279"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ORMDL1 number 17289"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OSBPL11 number 17295"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OSCP1 number 17304"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OSGEPL1 number 17306"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OSGIN2 number 17308"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OSR2 number 17311"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene OSTF1 number 17314"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OTC number 17317"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OTOGL number 17321"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OTOL1 number 17322"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OTOP3 number 17325"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OTUD3 number 17330"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OTUD5 number 17332"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OTUD6B number 17333"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OTULIN number 17336"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OVCA2 number 17340"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene OVCH2 number 17342"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene OVOL2 number 17344"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene P2RX3 number 17355"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene P2RY12 number 17360"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene P2RY13 number 17361"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene P2RY2 number 17363"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene P2RY6 number 17365"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene P2RY8 number 17366"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene P3H3 number 17369"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene P4HA3 number 17373"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene P4HB number 17374"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PAAF1 number 17377"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PABPC4 number 17380"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PAF1 number 17390"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PAH number 17397"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PAK1IP1 number 17402"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PAK4 number 17405"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PAM number 17414"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PAM16 number 17415"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PANK1 number 17419"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PANX1 number 17423"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PAQR4 number 17435"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PAQR5 number 17436"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PARM1 number 17449"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PARP1 number 17451"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PARP10 number 17452"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PARP14 number 17454"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PARP3 number 17457"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PARVA number 17463"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PARVG number 17465"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PAWR number 17472"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PAX9 number 17480"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PAXIP1 number 17482"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PBK number 17485"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PBRM1 number 17486"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PCBP3 number 17496"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PCCB number 17499"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PCDH8 number 17510"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PCID2 number 17519"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PCIF1 number 17520"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PCK1 number 17521"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PCLAF number 17523"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PCMTD2 number 17528"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PCNX4 number 17534"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PCOLCE2 number 17536"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PCP4 number 17538"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PCP4L1 number 17539"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PCSK1 number 17540"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PCSK9 number 17547"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PCYOX1L number 17549"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PCYT2 number 17552"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PDCD1 number 17555"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PDCD7 number 17563"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PDCL2 number 17565"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PDE6H number 17584"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PDE9A number 17589"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PDGFD number 17593"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PDHA1 number 17597"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PDIA2 number 17600"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PDILT number 17606"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PDK1 number 17607"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PDK4 number 17610"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PDLIM5 number 17615"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PDP2 number 17618"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PDX1 number 17627"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PDXK number 17629"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PDYN number 17630"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PECAM1 number 17647"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PECR number 17648"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PFDN1 number 17679"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PFDN5 number 17682"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PFKM number 17689"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PFN1 number 17690"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PGAP6 number 17699"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PGD number 17701"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PGGHG number 17703"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PGLS number 17706"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PGM2 number 17709"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PGM2L1 number 17710"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PHC3 number 17728"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PHF20 number 17739"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PHF8 number 17747"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PHKG2 number 17753"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PI3 number 17776"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PI4KB number 17780"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PIAS2 number 17782"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PIAS4 number 17784"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PIGG number 17799"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PIGN number 17804"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PIGR number 17807"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PIGT number 17809"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PIGY number 17814"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PIK3C2A number 17819"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PIK3R2 number 17829"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PIP4K2C number 17844"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PIPOX number 17851"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PIR number 17852"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PITPNB number 17857"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PJA2 number 17869"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PKDREJ number 17878"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PKIA number 17880"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PKLR number 17882"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PKP3 number 17892"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLA1A number 17894"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLA2G6 number 17902"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLAAT1 number 17906"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLBD1 number 17915"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PLCB2 number 17918"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLCB3 number 17919"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLCD4 number 17923"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLCG2 number 17926"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLCH1 number 17927"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLEK2 number 17942"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLEKHA2 number 17944"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLEKHB1 number 17951"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLEKHG5 number 17961"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLEKHH1 number 17963"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLEKHO2 number 17971"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLIN2 number 17975"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLK1 number 17976"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLPP5 number 17990"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLPP6 number 17991"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PLPPR3 number 17995"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PLXDC1 number 18002"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PMM2 number 18022"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PMS1 number 18026"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PNISR number 18029"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PNPLA6 number 18038"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PNRC2 number 18044"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene POC1B number 18046"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PODN number 18048"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene POGLUT1 number 18055"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene POGLUT2 number 18056"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene POLD3 number 18064"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene POLDIP3 number 18067"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene POLE number 18068"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene POLK number 18076"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene POLR1A number 18081"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene POLR3F number 18105"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene POM121 number 18111"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene POP1 number 18118"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene POU1F1 number 18127"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene POU3F3 number 18134"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene POU4F1 number 18136"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene POU5F1 number 18139"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPA1 number 18143"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPARGC1B number 18149"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPIL1 number 18172"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPIL6 number 18176"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPL number 18179"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPM1D number 18182"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPM1K number 18188"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPM1M number 18190"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPM1N number 18191"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PPOX number 18193"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPP1CC number 18196"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPP1R13B number 18202"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPP1R14A number 18204"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPP1R18 number 18213"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPP1R1B number 18215"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPP1R21 number 18218"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPP1R35 number 18222"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPP1R3A number 18225"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPP1R3C number 18227"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PPP1R8 number 18234"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPP1R9B number 18236"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PPP2R2D number 18242"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPP2R3A number 18243"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PPP2R5A number 18246"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPP2R5B number 18247"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPP3CB number 18252"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPP4R4 number 18260"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PPP6R1 number 18262"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PRAM1 number 18273"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PRDM2 number 18285"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PRDM4 number 18286"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PRDM9 number 18290"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PRDX4 number 18294"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PREB number 18297"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PREPL number 18302"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PREX1 number 18303"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PRG4 number 18305"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PRICKLE1 number 18306"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PRIMA1 number 18311"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PRIMPOL number 18312"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PRKAB2 number 18316"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PRKAG3 number 18320"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PRKAR2A number 18323"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PRKD3 number 18337"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PROX1 number 18369"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PRPS1 number 18387"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PRRC1 number 18404"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PRRT1 number 18412"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618392.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PRRX1 number 18417"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PRRX2 number 18418"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PRX number 18430"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PSD2 number 18436"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PSEN1 number 18439"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PSEN2 number 18440"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PSENEN number 18441"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PSMA4 number 18447"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PSMB1 number 18452"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PSMB2 number 18453"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PSMB6 number 18457"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PSMB7 number 18458"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PSMD12 number 18470"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PSMD13 number 18471"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PSMD5 number 18476"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PSMD9 number 18480"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PSPH number 18492"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTBP2 number 18501"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTCD3 number 18505"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTDSS1 number 18511"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTGER2 number 18519"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTGES3L number 18525"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTGFR number 18526"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTGIR number 18528"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTGS2 number 18533"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTK7 number 18541"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTN number 18544"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTP4A1 number 18545"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTPDC1 number 18549"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTPN11 number 18552"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTPN12 number 18553"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTPN13 number 18554"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTPN20 number 18558"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTPRN number 18577"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTPRS number 18582"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PTRH1 number 18586"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PUF60 number 18593"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PUM1 number 18594"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PURB number 18598"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PURG number 18599"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene PUS10 number 18600"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PUS7L number 18603"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PUSL1 number 18604"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PWP1 number 18607"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PWP2 number 18608"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PXDC1 number 18612"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PXN number 18618"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PYROXD1 number 18628"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene PYROXD2 number 18629"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene QARS1 number 18631"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene QPCT number 18634"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene QPRT number 18636"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene QSER1 number 18641"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene QSOX1 number 18642"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene R3HDML number 18651"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RAB12 number 18660"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RAB13 number 18661"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RAB1A number 18667"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RAB28 number 18678"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RAB38 number 18691"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RAB39A number 18692"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RAB3D number 18697"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RAB4A number 18707"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RAB4B number 18708"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RAB6A number 18713"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RAB7B number 18716"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RABEP1 number 18721"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RABEPK number 18723"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RAC2 number 18734"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RAD23B number 18744"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RADX number 18758"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RAG2 number 18762"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene RALBP1 number 18768"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RAMP2 number 18779"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RAP1GDS1 number 18792"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RAP2B number 18794"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene RASEF number 18821"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RASGEF1C number 18824"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RASGRP2 number 18828"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RASGRP3 number 18829"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RASL10B number 18832"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RASSF1 number 18836"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RB1CC1 number 18849"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RBBP7 number 18853"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RBBP8NL number 18855"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RBKS number 18863"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RBM10 number 18866"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RBM20 number 18875"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RBM25 number 18879"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RBM33 number 18883"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RBP5 number 18908"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene RBPJ number 18910"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RC3H1 number 18916"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RC3H2 number 18917"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RCAN1 number 18918"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RCAN3 number 18920"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RCC1 number 18922"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RCC1L number 18923"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RCC2 number 18924"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RCSD1 number 18934"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RD3L number 18937"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene RDH14 number 18939"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RDM1 number 18942"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RDX number 18943"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RECQL5 number 18948"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene REELD1 number 18949"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene REEP2 number 18951"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene REEP4 number 18953"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene REG4 number 18956"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RELCH number 18960"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RELT number 18964"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RERE number 18972"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RESF1 number 18975"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene RETSAT number 18981"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RFNG number 18999"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RFX4 number 19007"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RGL1 number 19014"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RGMB number 19018"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RGN number 19019"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RGS11 number 19024"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RGS13 number 19026"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RGS17 number 19029"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RGS19 number 19031"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RGS8 number 19040"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RHAG number 19043"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RHBDL3 number 19051"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RHOB number 19059"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene RHOT2 number 19070"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RIC8A number 19079"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RIMKLB number 19088"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RIN2 number 19094"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RINL number 19097"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RIOK1 number 19099"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RIOK2 number 19100"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RIOX2 number 19103"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RIPK3 number 19105"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RIT2 number 19114"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RMDN3 number 19122"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RMND1 number 19125"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RNASEH1 number 19128"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RNASET2 number 19132"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RND1 number 19133"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RNF111 number 19138"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RNF115 number 19140"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RNF121 number 19141"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RNF123 number 19143"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RNF125 number 19144"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RNF144B number 19153"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RNF152 number 19158"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RNF166 number 19161"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RNF175 number 19167"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RNF185 number 19171"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RNF217 number 19184"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RNF225 number 19188"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene RNF40 number 19198"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RNF7 number 19203"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ROMO1 number 19222"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene RP9 number 19233"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RPA2 number 19235"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene RPE65 number 19242"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RPF2 number 19244"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RPL13A number 19254"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene RPL35 number 19275"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RPL4 number 19283"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene RPL7 number 19286"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RPN1 number 19294"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RPN2 number 19295"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RPP25L number 19299"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RPP30 number 19300"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RPRM number 19306"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene RPRML number 19307"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene RPS17 number 19316"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RPS25 number 19324"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RPS27A number 19327"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene RPS27L number 19328"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RPS7 number 19346"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RPUSD4 number 19354"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RRAGD number 19357"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RRAS number 19358"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RRM2B number 19365"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RRP1B number 19370"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RRP8 number 19373"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RSPH9 number 19388"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RSPO2 number 19390"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RTCB number 19400"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RTN4IP1 number 19410"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RTN4R number 19411"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RTN4RL1 number 19412"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RUSF1 number 19431"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RWDD4 number 19437"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene RXFP2 number 19439"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene S100A1 number 19450"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene S100A14 number 19454"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene S100A6 number 19456"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene S100B number 19457"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene S1PR2 number 19460"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SAC3D1 number 19465"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SAFB number 19469"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SALL1 number 19471"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SAMD4B number 19484"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SAP30L number 19494"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SAPCD1 number 19495"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SAPCD2 number 19496"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SART1 number 19504"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SAV1 number 19513"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SAXO1 number 19514"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SBNO1 number 19522"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SBSPON number 19524"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SCAF11 number 19527"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SCAI number 19530"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SCAMP1 number 19531"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SCAMP4 number 19533"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SCARB1 number 19538"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SCN1B number 19556"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SCN2A number 19557"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SCRN3 number 19575"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SCUBE2 number 19581"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SCYL1 number 19584"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SDAD1 number 19587"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SDC1 number 19588"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SDHA number 19598"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SEC14L2 number 19615"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SEC22A number 19618"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SEC61G number 19632"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SEL1L number 19637"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SEL1L2 number 19638"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SELPLG number 19652"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SEMA4A number 19660"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SENP3 number 19675"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618375.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SEPTIN10 number 19683"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SEPTIN3 number 19687"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SEPTIN6 number 19689"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SEPTIN8 number 19691"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SERAC1 number 19693"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SERPINB5 number 19708"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SERPINE2 number 19711"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SERTAD1 number 19719"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SERTAD2 number 19720"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SERTM1 number 19723"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SESTD1 number 19728"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SETDB1 number 19740"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SEZ6 number 19742"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SF3A3 number 19748"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SF3B1 number 19749"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SFN number 19758"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SFRP1 number 19761"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SFT2D2 number 19765"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SFTPB number 19767"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SFTPD number 19768"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SGTA number 19795"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SH3GL3 number 19821"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SH3GLB2 number 19823"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SHFL number 19848"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SHH number 19849"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SHISA7 number 19855"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SHKBP1 number 19860"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SHOC1 number 19866"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SHROOM3 number 19875"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SHROOM4 number 19876"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SIAH3 number 19882"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SIGIRR number 19885"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SIM1 number 19894"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SIM2 number 19895"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SIN3A number 19897"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SIPA1 number 19900"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SIPA1L3 number 19903"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SIRT1 number 19904"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SIRT2 number 19905"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SIRT4 number 19907"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SIRT7 number 19910"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SKA2 number 19920"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SKIDA1 number 19925"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SKP2 number 19930"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLA2 number 19932"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLAMF1 number 19935"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLAMF8 number 19936"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SLC10A2 number 19938"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC10A3 number 19939"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SLC12A1 number 19944"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC12A8 number 19951"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC12A9 number 19952"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618375.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SLC13A4 number 19956"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC15A4 number 19961"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC16A9 number 19975"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC17A8 number 19979"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC18A1 number 19981"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC1A5 number 19990"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC22A4 number 20001"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC25A1 number 20009"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC25A10 number 20010"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC25A12 number 20012"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC25A14 number 20014"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC25A17 number 20017"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC25A20 number 20020"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC25A3 number 20030"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC25A30 number 20031"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC25A4 number 20040"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC25A40 number 20041"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC25A48 number 20048"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC26A4 number 20057"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC27A6 number 20067"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC29A3 number 20071"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC2A1 number 20073"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC30A3 number 20086"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC31A2 number 20094"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC35A1 number 20098"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC35A2 number 20099"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC35A3 number 20100"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC35B3 number 20105"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC35B4 number 20106"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC35C2 number 20108"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC35G1 number 20122"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC35G2 number 20123"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC37A2 number 20126"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC37A4 number 20128"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC38A2 number 20132"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC38A4 number 20134"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC38A7 number 20137"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC39A10 number 20141"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC39A13 number 20144"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC39A8 number 20152"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC39A9 number 20153"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC3A1 number 20154"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC49A3 number 20175"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC4A4 number 20182"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC4A5 number 20183"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC4A9 number 20186"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC5A11 number 20194"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC5A6 number 20199"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC6A15 number 20210"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC6A17 number 20212"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC6A20 number 20215"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC6A4 number 20216"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC6A5 number 20217"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC6A6 number 20218"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC7A1 number 20221"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC9A3 number 20241"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLC9A3R2 number 20243"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLCO1B3 number 20251"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SLF1 number 20258"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SLFNL1 number 20260"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SLX4 number 20275"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SMAD4 number 20279"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SMAD5 number 20280"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SMC6 number 20305"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SMCR8 number 20309"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SMG7 number 20314"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SMIM12 number 20319"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SMIM24 number 20327"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SMIM7 number 20338"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SMIM8 number 20339"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SMLR1 number 20341"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SMPDL3B number 20352"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SMPX number 20353"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SMURF1 number 20360"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SMYD2 number 20363"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SNAP25 number 20371"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SNCA number 20381"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SNCB number 20382"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SNPH number 20390"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SNRPB number 20401"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SNRPD1 number 20404"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SNRPF number 20408"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SNX11 number 20421"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SNX13 number 20423"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SNX16 number 20426"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SNX17 number 20427"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SNX2 number 20430"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SNX21 number 20432"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SNX3 number 20438"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SNX32 number 20441"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SNX6 number 20445"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SOCS1 number 20451"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SOD3 number 20460"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SOGA3 number 20462"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SORBS3 number 20466"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SOS2 number 20474"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SOWAHC number 20479"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SOX1 number 20481"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SOX12 number 20484"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SOX13 number 20485"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SOX30 number 20491"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SOX6 number 20494"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SP9 number 20506"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SPART number 20520"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPATA19 number 20527"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPATA24 number 20531"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPATA25 number 20532"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SPATA46 number 20535"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SPATA5 number 20537"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPATA7 number 20541"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPATS2L number 20546"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPCS3 number 20551"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPDL1 number 20553"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPECC1L number 20557"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPEG number 20560"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPIC number 20569"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPICE1 number 20570"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPINK4 number 20576"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPINT1 number 20577"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPINT2 number 20578"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPIRE2 number 20581"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPNS2 number 20583"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPO11 number 20585"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPP2 number 20595"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPR number 20599"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPRED2 number 20601"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPRED3 number 20602"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPRY2 number 20607"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPSB3 number 20613"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SPTSSB number 20625"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SRC number 20633"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SRGN number 20647"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SRL number 20649"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SRP72 number 20656"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SRRT number 20669"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618375.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SRSF10 number 20671"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SRSF11 number 20672"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SRSF2 number 20673"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SSB number 20683"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SSBP3 number 20686"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SSPN number 20694"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SST number 20701"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SSTR4 number 20705"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SSTR5 number 20706"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SSUH2 number 20708"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ST14 number 20710"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ST18 number 20711"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ST6GALNAC2 number 20721"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ST7 number 20726"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ST7L number 20727"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene STAB2 number 20735"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STAM2 number 20743"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STARD7 number 20754"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STAT3 number 20759"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STAT4 number 20760"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STAT6 number 20761"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STAU1 number 20762"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STIM2 number 20773"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STIP1 number 20776"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STK16 number 20779"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STK25 number 20783"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STK32B number 20788"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STK33 number 20790"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STK35 number 20791"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STK38 number 20793"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STKLD1 number 20798"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STOML2 number 20808"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STON1 number 20810"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STPG1 number 20814"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STPG3 number 20816"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene STPG4 number 20817"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STRA6 number 20818"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STRA8 number 20819"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STRADB number 20821"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STRC number 20822"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STRN4 number 20828"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STUM number 20833"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STX12 number 20836"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STX1A number 20841"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STX1B number 20842"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STX7 number 20848"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene STXBP4 number 20852"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SUCLG1 number 20862"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SUCO number 20865"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SULF1 number 20872"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SULT4A1 number 20876"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SUPT5H number 20888"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SUPT7L number 20890"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SUSD3 number 20896"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SV2B number 20904"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SVOPL number 20908"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SYBU number 20914"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SYDE2 number 20923"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SYMPK number 20926"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SYNC number 20929"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SYNDIG1 number 20930"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SYNGR1 number 20936"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SYNGR2 number 20937"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SYNGR4 number 20939"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene SYNJ1 number 20940"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SYT11 number 20953"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SYT8 number 20965"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SYTL1 number 20967"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SYTL2 number 20968"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene SZRD1 number 20973"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TAAR5 number 20976"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TADA1 number 20987"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TADA3 number 20990"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TAF10 number 20992"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TAF13 number 20995"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TAF4B number 21004"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TAF5 number 21005"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TAF8 number 21009"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TAL1 number 21018"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TAMALIN number 21021"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TANGO6 number 21026"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TAS1R1 number 21039"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TBC1D24 number 21069"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TBC1D30 number 21072"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TBC1D8 number 21077"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TBX19 number 21103"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TBX2 number 21104"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TBX20 number 21105"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TBX21 number 21106"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TBX3 number 21108"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TCEA2 number 21118"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TCEANC number 21120"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TCERG1 number 21122"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TCF15 number 21125"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TCF19 number 21126"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TCF20 number 21127"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TCF23 number 21129"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TCF7L1 number 21135"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TCIM number 21138"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TCOF1 number 21141"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TCTE1 number 21147"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TDP2 number 21154"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TDRD12 number 21157"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TDRD5 number 21160"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TDRKH number 21164"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TECPR1 number 21169"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TECPR2 number 21170"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TECR number 21171"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TEDC2 number 21175"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TES number 21204"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TESK1 number 21206"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TESMIN number 21208"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TESPA1 number 21209"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TEX15 number 21216"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TEX2 number 21217"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TEX30 number 21221"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TEX9 number 21231"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TFAP2A number 21233"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TFCP2L1 number 21242"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TFDP1 number 21243"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TFG number 21248"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TFRC number 21254"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TGDS number 21256"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TGFA number 21257"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TGFBI number 21262"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene THAP11 number 21278"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene THAP12 number 21279"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene THBS1 number 21289"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene THNSL1 number 21299"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene THOC2 number 21302"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene THY1 number 21319"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TIE1 number 21328"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TIGD5 number 21333"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TIGIT number 21334"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TIMM10 number 21336"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TIMM23B number 21343"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TIMM8A number 21347"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TKFC number 21364"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TKT number 21365"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TLCD3A number 21368"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TLCD4 number 21370"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TLE2 number 21373"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TLN2 number 21381"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TLNRD1 number 21382"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TLR5 number 21386"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TLR8 number 21388"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TLX1 number 21389"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TM4SF19 number 21397"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TM4SF4 number 21398"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TM6SF2 number 21401"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TM9SF1 number 21404"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TMBIM6 number 21411"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMCO4 number 21424"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMCO6 number 21425"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMED10 number 21427"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM102 number 21440"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618375.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TMEM120A number 21454"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM127 number 21460"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM129 number 21462"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM132E number 21469"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM139 number 21472"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TMEM164 number 21493"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TMEM167A number 21495"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM169 number 21498"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM17 number 21499"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM174 number 21503"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TMEM177 number 21505"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM184A number 21513"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM185A number 21516"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM186 number 21517"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TMEM192 number 21519"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM196 number 21520"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM198 number 21521"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM200A number 21523"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM203 number 21527"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TMEM204 number 21528"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM220 number 21542"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM230 number 21547"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM241 number 21556"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM254 number 21567"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM26 number 21572"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM260 number 21573"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM263 number 21575"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM41B number 21597"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM42 number 21598"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM45A number 21601"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM63B number 21616"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM67 number 21620"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM69 number 21621"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM71 number 21623"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM86A number 21631"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMEM8B number 21636"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMIE number 21645"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMPPE number 21653"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMPRSS13 number 21655"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMPRSS2 number 21657"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMPRSS4 number 21659"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TMTC3 number 21668"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMTC4 number 21669"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMX3 number 21674"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TMX4 number 21675"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TNC number 21676"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TNF number 21677"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TNFRSF13C number 21690"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TNFRSF1A number 21694"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TNFRSF1B number 21695"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TNFRSF21 number 21696"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TNFRSF25 number 21697"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TNFSF8 number 21708"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TNIP3 number 21713"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TNMD number 21718"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TNNC2 number 21721"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TNNT1 number 21725"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TNNT2 number 21726"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TNRC18 number 21732"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TOM1 number 21745"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TOM1L2 number 21747"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TOMM20 number 21748"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TOMM20L number 21749"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TOMM40 number 21752"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TONSL number 21758"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TOP2B number 21762"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TOP3B number 21764"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TOPAZ1 number 21765"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TOPBP1 number 21766"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TOPORS number 21767"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TOR2A number 21769"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TOR3A number 21770"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TOX4 number 21775"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TPH1 number 21795"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TPI1 number 21797"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TPM3 number 21801"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TPR number 21808"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TPRG1L number 21811"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TPST1 number 21814"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRADD number 21822"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRAF3 number 21825"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRAF3IP1 number 21826"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRAF3IP2 number 21827"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRAF4 number 21829"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRAM1 number 21836"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRAPPC10 number 21841"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRAPPC2 number 21846"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRAPPC3L number 21849"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRAPPC6B number 21853"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRARG1 number 21856"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TREX2 number 21862"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRIM2 number 21874"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRIM33 number 21881"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRIM45 number 21888"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRIM47 number 21890"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRIM56 number 21894"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRIM66 number 21899"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRIM69 number 21901"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRIM7 number 21902"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRIM72 number 21904"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRIP4 number 21912"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRIP6 number 21913"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRIR number 21915"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene TRMT61A number 21931"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRMT9B number 21933"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRMU number 21934"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRNAA-CGC number 21944"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAC-GCA_19 number 21993"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAD-GUC_9 number 22018"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAF-GAA_20 number 22064"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAF-GAA_31 number 22076"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618384.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAI-AAU_10 number 22152"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAI-UAU_13 number 22174"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618384.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAK-CUU number 22183"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAK-CUU_10 number 22185"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAL-AAG_3 number 22253"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAL-CAG_1 number 22267"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAM-CAU_32 number 22315"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618379.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAN-GUU_22 number 22338"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618379.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAN-GUU_23 number 22339"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618384.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAN-GUU_6 number 22343"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAR-CCG number 22411"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAR-CCU_5 number 22420"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618384.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAS-AGA_25 number 22451"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAS-CGA_15 number 22466"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618384.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAS-UGA_13 number 22501"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAS-UGA_16 number 22504"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAS-UGA_20 number 22509"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618384.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAT-UGU_3 number 22543"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAV-AAC_12 number 22554"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAV-AAC_23 number 22566"
<simpleError in plot.window(...): need finite 'ylim' values>
[1] "ERROR: Problem with gene TRNAV-AAC_3 number 22568"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRNAV-CAC_10 number 22577"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TRPM1 number 22644"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRPV2 number 22655"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TRPV4 number 22657"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TSC1 number 22661"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TSHB number 22674"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene TSHZ1 number 22676"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TSPAN15 number 22689"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TSPAN18 number 22691"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TSPAN19 number 22692"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TSPAN2 number 22693"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TSPAN3 number 22694"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TSPAN8 number 22702"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TSTD1 number 22715"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TSTD2 number 22716"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TTC13 number 22722"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TTC26 number 22732"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TTC39A number 22744"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TTL number 22761"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TTLL1 number 22762"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TTLL10 number 22763"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TTLL2 number 22766"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TTLL3 number 22767"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TTPAL number 22776"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TUBB4A number 22785"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TUBGCP4 number 22793"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TUFT1 number 22797"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TULP2 number 22799"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TXLNB number 22815"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TXLNG number 22816"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TXNDC15 number 22821"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TXNL4B number 22828"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TYRO3 number 22835"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TYRP1 number 22836"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene TYSND1 number 22837"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene U2AF1 number 22840"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UBA1 number 22846"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UBALD2 number 22856"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UBAP2 number 22859"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UBAP2L number 22860"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UBASH3B number 22862"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UBC number 22864"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene UBE2J1 number 22878"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UBE2L6 number 22882"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene UBE2R2 number 22889"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UBE2S number 22890"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UBE2U number 22892"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UBE3B number 22897"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UBIAD1 number 22903"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UBP1 number 22912"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UBR4 number 22918"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UCKL1 number 22937"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UFM1 number 22946"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UFSP1 number 22947"
<simpleError: scanVcf: scanVcf: scanTabix: 'NW_023618375.1' not present in tabix index
 path: inputs/Somatypus_Indels_final.vcf.gz
 index: inputs/Somatypus_Indels_final.vcf.gz.tbi
  path: inputs/Somatypus_Indels_final.vcf.gz>
[1] "Nothing to close."
[1] "ERROR: Problem with gene UGT8 number 22954"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UHMK1 number 22955"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UNC119B number 22969"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UNC13D number 22972"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UNC93A number 22982"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UNKL number 22987"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UPP1 number 22997"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene URB1 number 23000"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene URI1 number 23002"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UROS number 23006"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene USB1 number 23007"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene USF1 number 23009"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene USF2 number 23010"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene USH1C number 23012"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene USHBP1 number 23015"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene USP16 number 23024"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene USP36 number 23041"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene USP37 number 23042"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene USP38 number 23043"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene USP44 number 23048"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene USP45 number 23049"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene USP47 number 23051"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene UTP14A number 23065"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene UTP15 number 23066"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VAPA number 23089"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VASH2 number 23094"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VASN number 23095"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VCPKMT number 23108"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VEGFB number 23114"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene VEZF1 number 23118"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VGLL2 number 23121"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VIPAS39 number 23128"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VKORC1L1 number 23132"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VMAC number 23135"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene VOPP1 number 23137"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VPS18 number 23144"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VPS37B number 23157"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VPS37D number 23159"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene VPS39 number 23160"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VPS41 number 23161"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VPS4B number 23164"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VPS54 number 23169"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VRK1 number 23173"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VRK2 number 23174"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VSIG10L number 23179"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VWA3B number 23201"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VWA5B2 number 23203"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VWA7 number 23204"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VWCE number 23208"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene VWDE number 23209"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WASHC4 number 23222"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WASHC5 number 23223"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WBP4 number 23230"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WDFY4 number 23234"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WDR1 number 23236"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WDR18 number 23241"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WDR19 number 23242"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WDR20 number 23243"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WDR31 number 23249"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WDR37 number 23253"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WDR61 number 23270"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WDR64 number 23272"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WDR72 number 23275"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WDR93 number 23291"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WFDC1 number 23297"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WIF1 number 23303"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WNK1 number 23311"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WNT11 number 23318"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WNT7B number 23329"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WNT8B number 23331"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WRAP73 number 23335"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WRNIP1 number 23337"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WSCD2 number 23341"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WWP2 number 23348"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene WWTR1 number 23349"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene XAF1 number 23351"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene XCR1 number 23353"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene XG number 23355"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene XIRP1 number 23356"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene XPNPEP1 number 23367"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene XPO1 number 23370"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene XPO6 number 23373"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene XRCC4 number 23380"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene XRCC5 number 23381"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene XRN2 number 23384"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene YBX3 number 23398"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene YEATS4 number 23401"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene YIPF5 number 23409"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene YPEL1 number 23418"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene YPEL5 number 23422"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene YWHAG number 23431"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene YWHAZ number 23434"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZBED6CL number 23443"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZBTB24 number 23456"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZBTB32 number 23460"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZC3H11A number 23489"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZC3H12B number 23491"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZC3H12D number 23493"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZC3H18 number 23497"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZCCHC17 number 23507"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZCCHC7 number 23511"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZDHHC19 number 23526"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZDHHC2 number 23527"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZDHHC20 number 23528"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZDHHC23 number 23531"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZFAND4 number 23547"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZFP91 number 23559"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZFR number 23563"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZFX number 23564"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZFYVE1 number 23565"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZFYVE16 number 23566"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZFYVE9 number 23572"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZIC1 number 23578"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene ZIC3 number 23580"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZMAT2 number 23583"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZMAT3 number 23584"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZMYND15 number 23596"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZMYND8 number 23598"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZNF219 number 23608"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZNF335 number 23619"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZNF362 number 23622"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZNF410 number 23633"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZNF414 number 23634"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZNF436 number 23636"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene ZNF451 number 23638"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZNF488 number 23642"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZNF507 number 23644"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZNF511 number 23645"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZNF512 number 23646"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZNF576 number 23657"
<simpleError in apply(vars.vaf[, hosts], 1, function(vaf) {    all(vaf < MIN.VAF)}): dim(X) must have a positive length>
[1] "Nothing to close."
[1] "ERROR: Problem with gene ZNF622 number 23664"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZNF711 number 23682"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZNF740 number 23683"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZNRF1 number 23703"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZNRF3 number 23705"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZPBP2 number 23706"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZPLD1 number 23707"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZRANB3 number 23711"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZSWIM3 number 23715"
<simpleError in xy.coords(x, y, xlabel, ylabel, log): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZSWIM4 number 23716"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZSWIM9 number 23721"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZUP1 number 23722"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
[1] "ERROR: Problem with gene ZYG11B number 23727"
<simpleError in xy.coords(x, y): 'x' and 'y' lengths differ>
> 
