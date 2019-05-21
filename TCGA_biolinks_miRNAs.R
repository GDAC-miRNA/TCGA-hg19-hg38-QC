# 04 March 2019, 11h40 Pacific

# Using TCGAbiolinks to get TCGA-BLCA legacy and harmonized miRNA-seq data from the GDC, then comparing legacy vs harmonized RPM distributions
# Contrast mir-21 with the two members of the mir-24 family: hsa-mir-24-1 and -2.
# http://bioinformaticsfmrp.github.io/TCGAbiolinks/tcgaBiolinks.html

# Gordon Robertson, BC Cancer Agency, grobertson@bcgsc.ca


library(TCGAbiolinks) # v2.10.4
library(dplyr)


################################################################################
#-- set a working folder
################################################################################
# setwd("...")


################################################################################
#-- check which harmonized data release we'll query
################################################################################
getGDCInfo()$data_release
# "Data Release 15.0 - February 20, 2019"



################################################################################
#-- 1. Harmonized stem-loop data for the TCGA-BLCA cohort
################################################################################
#-- Query
query.mirna <- GDCquery(project = "TCGA-BLCA", 
                        data.category = "Transcriptome Profiling", 
                        data.type = "miRNA Expression Quantification",
                        workflow.type = "BCGSC miRNA Profiling",
                        sample.type = c("Primary solid Tumor"),
                        experimental.strategy = "miRNA-Seq"
)
#-- download
GDCdownload(query.mirna, 
            method = "api", 
            directory = "GDCdata_hg38_stemloops", 
            files.per.chunk = 50 )
# Downloading data for project TCGA-BLCA
# GDCdownload will download 417 files. A total of 20.98574 MB
# Downloading chunk 1 of 9 (50 files, size = 2.516258 MB) as Sun_Feb_24_09_50_40_2019_0.tar.gz

#-- Prepare: save an .RData file into working folder, then delete all intermediate downloadeded files
# When using the funciton GDCprepare there is an argument called SummarizedExperiment which defines the output type a Summarized Experiment (default option) or a data frame. 
# To create a summarized Experiment object we annotate the data with genomic positions with last patch release version of the genome available. 
# For legacy data (data aligned to hg19) TCGAbiolinks is using GRCh37.p13 and for harmonized data (data aligned to hg38) now it is using GRCh38.p7 (May 2017).
hg38_stemloop_data <- GDCprepare(query.mirna, 
                                 save = T, 
                                 save.filename = "hg38_mirna_stemloops.RData", 
                                 directory = "GDCdata_hg38_stemloops", 
                                 remove.files.prepared = TRUE 
                                 )


#-- Now load the .RData file
hg38_stemloop_datafile <- get(load("hg38_mirna_stemloops.RData"))

class(hg38_stemloop_datafile)
# "data.frame"
dim(hg38_stemloop_datafile)
# 1881 1252

#-- For each sample, the file has a read_count, RPM, and crossmapped column
hg38_stemloop_datafile[1:5,1:5]

#-- get only the miRNA stemloop name and the RPM columns
hg38_stemloop_RPMs <- hg38_stemloop_datafile[ , c(1,grep("million", colnames(hg38_stemloop_datafile)))]
dim(hg38_stemloop_RPMs)
# 1881 418

#-- check colnames
table(substr(colnames(hg38_stemloop_RPMs), 1,18))
# miRNA_ID reads_per_million_ 
#        1                417 

#-- shorten colnames
colnames(hg38_stemloop_RPMs) <- sub("reads_per_million_miRNA_mapped_", "", colnames(hg38_stemloop_RPMs))
hg38_stemloop_RPMs[1:5,1:5]
#       miRNA_ID TCGA-FD-A6TE-01A-12R-A33A-13 TCGA-C4-A0F7-01A-11R-A085-13 TCGA-GU-AATO-01A-11R-A39B-13 TCGA-K4-A83P-01A-11R-A358-13
# 1 hsa-let-7a-1                    9736.1845                     4526.388                    6204.6717                     5872.616
# 2 hsa-let-7a-2                    9634.7070                     4398.605                    6190.3701                     5843.260
# 3 hsa-let-7a-3                    9827.2697                     4500.592                    6289.6559                     5788.741
# 4   hsa-let-7b                   11958.4499                    11669.622                    7707.9844                     5652.795
# 5   hsa-let-7c                     205.8587                     1211.836                     908.9734                     1540.671


#-- get only primary tumours
hg38_primary_tumour_barcodes <- colnames(hg38_stemloop_RPMs)[substr(colnames(hg38_stemloop_RPMs), 14,15) == "01"]
length(hg38_primary_tumour_barcodes)
# 417
head(hg38_primary_tumour_barcodes)

#-- get the columns
hg38_stemloop_RPMs_01 <- hg38_stemloop_RPMs[ , c("miRNA_ID", hg38_primary_tumour_barcodes)]
dim(hg38_stemloop_RPMs_01)
# 1881 418

#-- Look at what we have
hg38_stemloop_RPMs_01[1:5,1:5]



################################################################################
#-- hg38 stemloops: get hsa-mir-24-1 and 24-2 RPMs
################################################################################
grep("-24-", hg38_stemloop_RPMs_01$miRNA_ID, value = T)
# "hsa-mir-24-1" "hsa-mir-24-2"

#-- hsa-mir-24-1
hg38_hsa.mir.24.1_stemloop_01_RPMs <- as.numeric(subset(hg38_stemloop_RPMs_01, miRNA_ID == "hsa-mir-24-1"))

#-- remove first element's NA
is.na(hg38_hsa.mir.24.1_stemloop_01_RPMs)
hg38_hsa.mir.24.1_stemloop_01_RPMs <- hg38_hsa.mir.24.1_stemloop_01_RPMs[-1]
sort(hg38_hsa.mir.24.1_stemloop_01_RPMs)
#-- quick log2 boxplot
boxplot(log2(hg38_hsa.mir.24.1_stemloop_01_RPMs), ylim = c(0,20))


#-- hsa-mir-24-2 stem-loop
hg38_hsa.mir.24.2_stemloop_01_RPMs <- as.numeric(subset(hg38_stemloop_RPMs, miRNA_ID == "hsa-mir-24-2"))
#-- remove first element NA
hg38_hsa.mir.24.2_stemloop_01_RPMs <- hg38_hsa.mir.24.2_stemloop_01_RPMs[-1]
head(hg38_hsa.mir.24.2_stemloop_01_RPMs)


#-- compare them
boxplot(log2(hg38_hsa.mir.24.1_stemloop_01_RPMs), log2(hg38_hsa.mir.24.2_stemloop_01_RPMs), log = "y", names = c("24-1","24-2"), las = 1)
wilcox.test(hg38_hsa.mir.24.1_stemloop_01_RPMs, hg38_hsa.mir.24.2_stemloop_01_RPMs)
# Wilcoxon rank sum test with continuity correction
# W = 86477, p-value = 0.8932



################################################################################
#-- The miR-21 stem-loop is MI0000077
################################################################################
grep("mir-21", hg38_stemloop_RPMs$miRNA_ID, value = T)
# "hsa-mir-21"

hg38_hsa.mir.21_RPMs_stemloop_01 <- as.numeric(subset(hg38_stemloop_RPMs, miRNA_ID == "hsa-mir-21"))
is.na(hg38_hsa.mir.21_RPMs_stemloop_01)

#-- remove first element NA
hg38_hsa.mir.21_RPMs_stemloop_01 <- hg38_hsa.mir.21_RPMs_stemloop_01[-1]
sort(hg38_hsa.mir.21_RPMs_stemloop_01)

#-- quick look
boxplot(log2(hg38_hsa.mir.21_RPMs_stemloop_01), ylim = c(0,20))





################################################################################
#-- 2. Harmonized hg38 isoform data
################################################################################
query.mirna.isoform <- GDCquery(project = "TCGA-BLCA", 
                        data.category = "Transcriptome Profiling", 
                        data.type = "Isoform Expression Quantification",
                        workflow.type = "BCGSC miRNA Profiling",
                        experimental.strategy = "miRNA-Seq",
                        sample.type = c("Primary solid Tumor")
                        )

GDCdownload(query.mirna.isoform, 
            method = "api", 
            directory = "GDCdata_hg38_isoforms", 
            files.per.chunk = 50)

#-- Prepare: save an .RData file into working folder, then delete all intermediate downloadeded files
hg38_isoform_data <- GDCprepare(query.mirna.isoform, 
                                directory = "GDCdata_hg38_isoforms", 
                                save = T, 
                                save.filename = "hg38_mirna_isoforms.RData",
                                remove.files.prepared = TRUE 
                                )

#-- load the .RData from the working directory
hg38_isoform_datafile <- get(load("hg38_mirna_isoforms.RData"))

class(hg38_isoform_datafile)
# "spec_tbl_df" "tbl_df"      "tbl"         "data.frame" 
dim(hg38_isoform_datafile)
# 2025288       7

# For each sample, the file has a read_count, RPM, and crossmapped column
head(hg38_isoform_datafile)
# A tibble: 6 x 7
#   miRNA_ID     isoform_coords                read_count reads_per_million_miRNA_mapped `cross-mapped` miRNA_region        barcode                     
#   <chr>        <chr>                              <int>                          <dbl> <chr>          <chr>               <chr>                       
# 1 hsa-let-7a-1 hg38:chr9:94175961-94175982:+          1                          0.354 N              mature,MIMAT0000062 TCGA-E7-A7PW-01A-11R-A358-13
# 2 hsa-let-7a-1 hg38:chr9:94175961-94175983:+          4                          1.41  N              mature,MIMAT0000062 TCGA-E7-A7PW-01A-11R-A358-13
# 3 hsa-let-7a-1 hg38:chr9:94175961-94175984:+         13                          4.60  N              mature,MIMAT0000062 TCGA-E7-A7PW-01A-11R-A358-13
# 4 hsa-let-7a-1 hg38:chr9:94175961-94175985:+          1                          0.354 N              mature,MIMAT0000062 TCGA-E7-A7PW-01A-11R-A358-13
# 5 hsa-let-7a-1 hg38:chr9:94175962-94175981:+         33                         11.7   N              mature,MIMAT0000062 TCGA-E7-A7PW-01A-11R-A358-13
# 6 hsa-let-7a-1 hg38:chr9:94175962-94175982:+       1662                        588.    N              mature,MIMAT0000062 TCGA-E7-A7PW-01A-11R-A358-13


################################################################################
#-- hsa-miR-24-3p is MIMAT0000080
################################################################################
hg38_MIMAT0000080 <- hg38_isoform_datafile[grep("MIMAT0000080", hg38_isoform_datafile$miRNA_region),]
dim(hg38_MIMAT0000080)
# 15600     7
head(hg38_MIMAT0000080)

#-- group by barcode, and get the sum of RPMs for all isoforms for a barcode
hg38_MIMAT0000080_01_grouped_sums <- hg38_MIMAT0000080 %>%
  filter(
    substring(barcode,14,15) == "01"
  ) %>%
  group_by(barcode) %>%
  summarize(sum = sum(reads_per_million_miRNA_mapped))
hg38_MIMAT0000080_01_grouped_sums

#-- check primary tumours
table(substr(hg38_MIMAT0000080_01_grouped_sums$barcode, 14,15))
# 01 417

range(hg38_MIMAT0000080_01_grouped_sums$sum)
# 408.6798 19696.0160

#-- quick look at an EDF with a log-x scale
plot(ecdf(hg38_MIMAT0000080_01_grouped_sums$sum), log = "x", xlim = c(400,20000))


################################################################################
#-- hsa-miR-24-1-5p is MIMAT0000079
################################################################################
hg38_MIMAT0000079 <- hg38_isoform_datafile[grep("MIMAT0000079", hg38_isoform_datafile$miRNA_region),]
dim(hg38_MIMAT0000079)
# 4182     7
head(hg38_MIMAT0000079)

hg38_MIMAT0000079_01_grouped_sums <- hg38_MIMAT0000079 %>%
  filter(
    substring(barcode,14,15) == "01"
  ) %>%
  group_by(barcode) %>%
  summarize(sum = sum(reads_per_million_miRNA_mapped))
hg38_MIMAT0000079_01_grouped_sums


################################################################################
#-- hsa-miR-24-2-5p is MIMAT0004497
################################################################################
hg38_MIMAT0004497 <- hg38_isoform_datafile[grep("MIMAT0004497", hg38_isoform_datafile$miRNA_region),]
dim(hg38_MIMAT0004497)
# 4081     7
head(hg38_MIMAT0004497)

hg38_MIMAT0004497_01_grouped_sums <- hg38_MIMAT0004497 %>%
  filter(
    substring(barcode,14,15) == "01"
  ) %>%
  group_by(barcode) %>%
  summarize(sum = sum(reads_per_million_miRNA_mapped))
hg38_MIMAT0004497_01_grouped_sums
#   barcode                        sum
#   <chr>                        <dbl>
# 1 TCGA-2F-A9KO-01A-11R-A38M-13  43.6
# 2 TCGA-2F-A9KP-01A-11R-A38M-13  19.9
# 3 TCGA-2F-A9KQ-01A-11R-A38M-13  22.4
# 4 TCGA-2F-A9KR-01A-11R-A38M-13  31.7


boxplot(
  log2(hg38_MIMAT0000079_01_grouped_sums$sum), # 24-1-5p
  log2(hg38_MIMAT0000080_01_grouped_sums$sum), # 24-3p
  log2(hg38_MIMAT0004497_01_grouped_sums$sum), # 24-2-5p
  ylim=c(0,20),
  col = c("orange", "red", "skyblue2"),
  notch = T,
  las = 1,
)





################################################################################
#-- 3. Legacy hg19 stem-loops, boxplots
################################################################################
query.legacy.mirna <- GDCquery(project = "TCGA-BLCA", 
                               legacy= TRUE,
                               data.category = "Gene expression", 
                               data.type = "miRNA gene quantification",
                               sample.type = c("Primary solid Tumor","Solid Tissue Normal"),
                               file.type = "hg19.mirbase20")
GDCdownload(query.legacy.mirna, 
            method = "api", 
            directory = "GDCdata_hg19_stemloops", 
            files.per.chunk = 50)

#-- Prepare: save an .RData file into working folder, then delete all intermediate downloadeded files
hg19_stemloop_data <- GDCprepare(query.legacy.mirna, 
                                 directory = "GDCdata_hg19_stemloops",
                                 save = T, 
                                 save.filename = "hg19_mirna_stemloops.RData",
                                 remove.files.prepared = TRUE 
                                 )

#-- load the legacy .RData from the working folder
hg19_stemloop_datafile <- get(load("hg19_mirna_stemloops.RData"))

class(hg19_stemloop_datafile)
# "data.frame"
dim(hg19_stemloop_datafile)
# 1870 1285
# For each sample, the file has a read_count, RPM, and crossmapped column
hg19_stemloop_datafile[1:5,1:5]

#-- get only the miRNA stemloop name and the RPM columns
hg19_stemloop_RPMs <- hg19_stemloop_datafile[ , c(1,grep("million", colnames(hg19_stemloop_datafile)))]
dim(hg19_stemloop_RPMs)
# 1870 429

#-- check colnames
table(substr(colnames(hg19_stemloop_RPMs), 1,18))
# miRNA_ID reads_per_million_ 
#        1                428

#-- shorten the colnames 
colnames(hg19_stemloop_RPMs) <- sub("reads_per_million_miRNA_mapped_", "", colnames(hg19_stemloop_RPMs))
hg19_stemloop_RPMs[1:5,1:5]
#       miRNA_ID TCGA-FD-A6TE-01A-12R-A33A-13 TCGA-C4-A0F7-01A-11R-A085-13 TCGA-GU-AATO-01A-11R-A39B-13 TCGA-K4-A83P-01A-11R-A358-13
# 1 hsa-let-7a-1                    9736.1845                     4526.388                    6204.6717                     5872.616
# 2 hsa-let-7a-2                    9634.7070                     4398.605                    6190.3701                     5843.260
# 3 hsa-let-7a-3                    9827.2697                     4500.592                    6289.6559                     5788.741
# 4   hsa-let-7b                   11958.4499                    11669.622                    7707.9844                     5652.795
# 5   hsa-let-7c                     205.8587                     1211.836                     908.9734                     1540.671

#-- check only primary tumours
table(substr(colnames(hg19_stemloop_RPMs), 14,15))
#    01  11 
# 1 409  19 

#-- get only primary tumours
hg19_primary_tumour_barcodes <- colnames(hg19_stemloop_RPMs)[substr(colnames(hg19_stemloop_RPMs), 14,15) == "01"]
length(hg19_primary_tumour_barcodes)
# 409
head(hg19_primary_tumour_barcodes)

hg19_stemloop_RPMs_01 <- hg19_stemloop_RPMs[ , c("miRNA_ID", hg19_primary_tumour_barcodes)]
dim(hg19_stemloop_RPMs_01)
# 1870 410
hg19_stemloop_RPMs_01[1:5,1:5]



################################################################################
#-- get hsa-mir-24-1 and 24-2 stem-loop RPMs
################################################################################
grep("mir-24", hg19_stemloop_RPMs_01$miRNA_ID, value = T)
# "hsa-mir-24-1" "hsa-mir-24-2" "hsa-mir-2467"

#-- mir-24-1 stem-loop
hg19_hsa.mir.24.1_stemloop_01_RPMs <- as.numeric(subset(hg19_stemloop_RPMs_01, miRNA_ID == "hsa-mir-24-1"))
is.na(hg19_hsa.mir.24.1_stemloop_01_RPMs)
hg19_hsa.mir.24.1_stemloop_01_RPMs[1:5]
#-- remove the initial NA
hg19_hsa.mir.24.1_stemloop_01_RPMs <- hg19_hsa.mir.24.1_stemloop_01_RPMs[-1]
hg19_hsa.mir.24.1_stemloop_01_RPMs[1:5]
hg19_hsa.mir.24.1_stemloop_01_RPMs <- hg19_hsa.mir.24.1_stemloop_01_RPMs[!is.na(hg19_hsa.mir.24.1_stemloop_01_RPMs)]
sort(hg19_hsa.mir.24.1_stemloop_01_RPMs)
boxplot(hg19_hsa.mir.24.1_stemloop_01_RPMs)

#-- mir-24-2 stem-loop
hg19_hsa.mir.24.2_stemloop_01_RPMs <- as.numeric(subset(hg19_stemloop_RPMs_01, miRNA_ID == "hsa-mir-24-2"))
hg19_hsa.mir.24.2_stemloop_01_RPMs <- hg19_hsa.mir.24.2_stemloop_01_RPMs[-1]

#hsa.mir.24.1_and_2_RPMs_hg19 <- c(hg19_hsa.mir.24.1_stemloop_01_RPMs, hg19_hsa.mir.24.2_stemloop_01_RPMs)
#sort(hsa.mir.24.1_and_2_RPMs_hg19)


#-- compare them, boxplot
#pdf("miR-24_stem-loops.legacy_vs_harmonized.pdf", height = 3., width = 3.25)
par(mar = c(3.25,4,1.5,0.5), mgp = c(3,0.4,0), tck = -0.02)
boxplot(hg19_hsa.mir.24.1_stemloop_01_RPMs, 
        hg19_hsa.mir.24.2_stemloop_01_RPMs, 
        hg38_hsa.mir.24.1_stemloop_01_RPMs, 
        hg38_hsa.mir.24.2_stemloop_01_RPMs, 
        log = "y", las = 1,
        names = c("24-1", "24-2", "24-1", "24-2"),
        col = c("grey85","grey85","skyblue2","skyblue2"),
        cex.axis = 0.9,
        notch = T
        )
grid(nx = NA, ny = NULL, lty = 1, lwd = 0.5)
mtext("GRCh37      GRCh38", side = 1, line = 1.75, cex = 1.2)
mtext("Stem-loop RPM", side = 2, line = 2.75, cex = 1.2)
mtext("hsa-mir-24-1/2, BLCA", side = 3, line = 0.2, cex = 1.2)
#dev.off()

wilcox.test(hg38_hsa.mir.24.1_stemloop_01_RPMs, hg38_hsa.mir.24.2_stemloop_01_RPMs)
# Wilcoxon rank sum test with continuity correction
# W = 86477, p-value = 0.8932




################################################################################
#-- miR-21 stem-loop is MI0000077
################################################################################
hg19_stemloop_RPMs_01[1:5,1:5]
table(substr(colnames(hg19_stemloop_RPMs_01),14,15))

grep("mir-21", hg19_stemloop_RPMs_01$miRNA_ID, value = T)
# "hsa-mir-21" ...

hg19_hsa.mir.21_stemloop_RPMs_01 <- as.numeric(subset(hg19_stemloop_RPMs_01, miRNA_ID == "hsa-mir-21"))
hg19_hsa.mir.21_stemloop_RPMs_01[1:5,1:5]
is.na(hg19_hsa.mir.21_stemloop_RPMs_01)

#-- remove first element NA
hg19_hsa.mir.21_stemloop_RPMs_01 <- hg19_hsa.mir.21_stemloop_RPMs_01[-1]
sort(hg19_hsa.mir.21_stemloop_RPMs_01)
boxplot(hg19_hsa.mir.21_stemloop_RPMs_01, hg38_hsa.mir.21_RPMs_stemloop_01)

#-- medians?
median(hg19_hsa.mir.21_stemloop_RPMs_01/1000)
# 301.45
median(hg38_hsa.mir.21_RPMs_stemloop_01/1000)
# 301.8099

# Wilcoxon rank sum test with continuity correction
wilcox.test(hg19_hsa.mir.21_stemloop_RPMs_01, hg38_hsa.mir.21_RPMs_stemloop_01)
# W = 85905, p-value = 0.8547
#alternative hypothesis: true location shift is not equal to 0

#-- Ranges of values
range(hsa.mir.21_RPMs_hg19/1000)
# 3.803309 639.675890
range(hsa.mir.21_stemloop_01_RPMs/1000)
# 27.87122 636.56583

#-- Boxplot
#pdf("hsa-miR-21_stemloop.legacy_vs_harmonized.pdf", height = 3., width = 2.75)
par(mar = c(2.,3.75,1.5,0.5), mgp = c(3,0.4,0), tck = -0.02)
boxplot(log2(hsa.mir.21_RPMs_hg19), 
        log2(hsa.mir.21_stemloop_01_RPMs), 
        log = "y", las = 1,
        names = c("GRCh37", "GRCh38"),
        col = c("grey85","skyblue2"),
        cex.axis = 0.9,
        notch = T,
        xaxt = "n", yaxt = "n"
)
#grid(nx = NA, ny = NULL, lty = 1, lwd = 0.5)
for(i in c(5,10,20,50,100,200,500,1000)) lines(c(0,3), c(i,i), lty = 1, lwd = 0.5, col = "grey")
axis(side = 1, tick = F, labels = c("GRCh37","GRCh38"), at = 1:2, line = 0.0, cex.axis = 1.1)
axis(side = 2, tick = F, labels = c(5,10,20,50,100,200,500,1000), at = c(5,10,20,50,100,200,500,1000), line = -0.2, cex.axis = 1., las = 1)
mtext("log2(RPM)", side = 2, line = 2.25, cex = 1.3)
mtext("hsa-mir-21, BLCA", side = 3, line = 0.2, cex = 1.2)
#dev.off()




################################################################################
#-- 4. Legacy hg19 isoform data
################################################################################
query.legacy.mirna.isoform <- GDCquery(project = "TCGA-BLCA",
                                       legacy= TRUE,
                                       data.category = "Gene expression", 
                                       data.type = "miRNA isoform quantification",
                                       sample.type = c("Primary solid Tumor","Solid Tissue Normal"),
                                       file.type = "hg19.isoform")

GDCdownload(query.legacy.mirna.isoform, 
            method = "api", 
            directory = "GDCdata_hg19_isoforms", 
            files.per.chunk = 50)

#-- Prepare: save an .RData file into working folder, then delete all intermediate downloadeded files
hg19_isoform_data <- GDCprepare(query.legacy.mirna.isoform, 
                                directory = "GDCdata_hg19_isoforms",
                                save = T, 
                                save.filename = "hg19_mirna_isoforms.RData",
                                remove.files.prepared = TRUE 
)


#-- load the .RData from the working folder
hg19_isoform_datafile <- get(load("hg19_mirna_isoforms.RData"))

class(hg19_isoform_datafile)
# "data.frame" 
dim(hg19_isoform_datafile)
# 2077720       7

# For each sample, the file has a read_count, RPM, and crossmapped column
head(hg19_isoform_datafile)
#       miRNA_ID             isoform_coords read_count reads_per_million_miRNA_mapped cross-mapped        miRNA_region                      barcode
# 1 hsa-let-7a-1 hg19:9:96938243-96938264:+          1                       0.142427            N mature,MIMAT0000062 TCGA-CU-A0YN-11A-11R-A10V-13
# 2 hsa-let-7a-1 hg19:9:96938243-96938265:+          1                       0.142427            N mature,MIMAT0000062 TCGA-CU-A0YN-11A-11R-A10V-13
# 3 hsa-let-7a-1 hg19:9:96938243-96938266:+          2                       0.284855            N mature,MIMAT0000062 TCGA-CU-A0YN-11A-11R-A10V-13
# 4 hsa-let-7a-1 hg19:9:96938244-96938263:+        105                      14.954877            N mature,MIMAT0000062 TCGA-CU-A0YN-11A-11R-A10V-13
# 5 hsa-let-7a-1 hg19:9:96938244-96938264:+       2801                     398.939144            N mature,MIMAT0000062 TCGA-CU-A0YN-11A-11R-A10V-13
# 6 hsa-let-7a-1 hg19:9:96938244-96938265:+       6412                     913.244480            N mature,MIMAT0000062 TCGA-CU-A0YN-11A-11R-A10V-13



################################################################################
#-- miR-24-1-5p 
################################################################################
hg19_MIMAT0000079 <- hg19_isoform_datafile[grep("MIMAT0000079", hg19_isoform_datafile$miRNA_region),]
dim(hg19_MIMAT0000079)
# 4365     7
head(hg19_MIMAT0000079)

#-- get only the primary tumour samples
hg19_MIMAT0000079_01_grouped_sums <- hg19_MIMAT0000079 %>%
  filter(
    substring(barcode,14,15) == "01"
  ) %>%
  group_by(barcode) %>%
  summarize(sum = sum(reads_per_million_miRNA_mapped))


################################################################################
#-- miR-24-2-5p
################################################################################
hg19_MIMAT0004497 <- hg19_isoform_datafile[grep("MIMAT0004497", hg19_isoform_datafile$miRNA_region),]
dim(hg19_MIMAT0004497)
# 4219     7
head(hg19_MIMAT0004497)

#-- get only the primary tumour samples
hg19_MIMAT0004497_01_grouped_sums <- hg19_MIMAT0004497 %>%
  filter(
    substring(barcode,14,15) == "01"
  ) %>%
  group_by(barcode) %>%
  summarize(sum = sum(reads_per_million_miRNA_mapped))


################################################################################
#-- miR-24-3p can be expressed from hsa-mir-24-1 and 24-2 pri-miRNAs
################################################################################
hg19_MIMAT0000080 <- hg19_isoform_datafile[grep("MIMAT0000080", hg19_isoform_datafile$miRNA_region),]
dim(hg19_MIMAT0000080)
# 9927     7
head(hg19_MIMAT0000080)

#-- get only the primary tumour samples
hg19_MIMAT0000080_01_grouped_sums <- hg19_MIMAT0000080 %>%
  filter(
    substring(barcode,14,15) == "01"
  ) %>%
  group_by(barcode) %>%
  summarize(sum = sum(reads_per_million_miRNA_mapped))

dim(hg19_MIMAT0000080_01_grouped_sums)
# 409     7

table(substr(hg19_MIMAT0000080_01_grouped_sums$barcode, 14,15))
# 01 = 409

range(hg19_MIMAT0000080_01_grouped_sums$sum)
# 407.1555 19614.1421

#-- plot an EDF with a log-x scale
plot(ecdf(hg19_MIMAT0000080_01_grouped_sums$sum), log = "x", xlim = c(400,20000))


#-- Wilcoxon rank sum test with continuity correction
wilcox.test(hg19_MIMAT0000080_01_grouped_sums$sum, 
            hg38_MIMAT0000080_01_grouped_sums$sum)
# W = 84644, p-value = 0.8537
# alternative hypothesis: true location shift is not equal to 0

#-- medians?
median(hg19_MIMAT0000080_01_grouped_sums$sum)
# 2199.587
median(hg38_MIMAT0000080_01_grouped_sums$sum)
# 2202.07
100*((median(hg38_MIMAT0000080_01_grouped_sums$sum) - median(hg19_MIMAT0000080_01_grouped_sums$sum))/median(hg19_MIMAT0000080_01_grouped_sums$sum))
# 0.1129





################################################################################
#-- hsa-miR-21-5p is MIMAT0000076; hsa-miR-21-3p is MIMAT0004494
################################################################################

#-- Harmonized hg38 - miR-21-5p
hg38_MIMAT0000076 <- hg38_isoform_datafile[grep("MIMAT0000076", hg38_isoform_datafile$miRNA_region),]
dim(hg38_MIMAT0000076)
# 21150     7
head(hg38_MIMAT0000076)
# miRNA_ID   isoform_coords                 read_count reads_per_million_miRNA_mapped `cross-mapped` miRNA_region        barcode                     
# <chr>      <chr>                               <int>                          <dbl> <chr>          <chr>               <chr>                       
# 1 hsa-mir-21 hg38:chr17:59841270-59841293:+          1                          0.354 N              mature,MIMAT0000076 TCGA-E7-A7PW-01A-11R-A358-13
# 2 hsa-mir-21 hg38:chr17:59841272-59841293:+          6                          2.12  N              mature,MIMAT0000076 TCGA-E7-A7PW-01A-11R-A358-13
# 3 hsa-mir-21 hg38:chr17:59841272-59841294:+        120                         42.4   N              mature,MIMAT0000076 TCGA-E7-A7PW-01A-11R-A358-13
# 4 hsa-mir-21 hg38:chr17:59841272-59841295:+        568                        201.    N              mature,MIMAT0000076 TCGA-E7-A7PW-01A-11R-A358-13
# 5 hsa-mir-21 hg38:chr17:59841272-59841296:+         13                          4.60  N              mature,MIMAT0000076 TCGA-E7-A7PW-01A-11R-A358-13
# 6 hsa-mir-21 hg38:chr17:59841272-59841297:+          2                          0.707 N              mature,MIMAT0000076 TCGA-E7-A7PW-01A-11R-A358-13

#-- get only the primary tumour samples
hg38_MIMAT0000076_01_grouped_sums <- hg38_MIMAT0000076 %>%
  filter(
    substring(barcode,14,15) == "01"
  ) %>%
  group_by(barcode) %>%
  summarize(sum = sum(reads_per_million_miRNA_mapped))

range(hg38_MIMAT0000076_01_grouped_sums$sum)
# 27548.98 632169.16

table(substr(hg38_MIMAT0000076_01_grouped_sums$barcode, 14,15))
# 01 417


################################################################################
#-- hg38 hsa-miR-21-3p is MIMAT0004494
################################################################################
hg38_MIMAT0004494 <- hg38_isoform_datafile[grep("MIMAT0004494", hg38_isoform_datafile$miRNA_region),]
dim(hg38_MIMAT0004494)
# 4081     7
head(hg38_MIMAT0004494)

hg38_MIMAT0004494_01_grouped_sums <- hg38_MIMAT0004494 %>%
  filter(
    substring(barcode,14,15) == "01"
  ) %>%
  group_by(barcode) %>%
  summarize(sum = sum(reads_per_million_miRNA_mapped))
hg38_MIMAT0004494_01_grouped_sums



################################################################################
#-- Legacy hg19 - hsa-miR-21-5p
################################################################################
hg19_MIMAT0000076 <- hg19_isoform_datafile[grep("MIMAT0000076", hg19_isoform_datafile$miRNA_region),]
dim(hg19_MIMAT0000076)
# 21541     7
head(hg19_MIMAT0000076)
# miRNA_ID              isoform_coords read_count reads_per_million_miRNA_mapped cross-mapped        miRNA_region                      barcode
# 1333 hsa-mir-21 hg19:17:57918633-57918652:+          1                       0.142427            N mature,MIMAT0000076 TCGA-CU-A0YN-11A-11R-A10V-13
# 1334 hsa-mir-21 hg19:17:57918633-57918653:+          1                       0.142427            N mature,MIMAT0000076 TCGA-CU-A0YN-11A-11R-A10V-13
# 1335 hsa-mir-21 hg19:17:57918633-57918654:+          1                       0.142427            N mature,MIMAT0000076 TCGA-CU-A0YN-11A-11R-A10V-13
# 1336 hsa-mir-21 hg19:17:57918633-57918655:+         60                       8.545644            N mature,MIMAT0000076 TCGA-CU-A0YN-11A-11R-A10V-13
# 1337 hsa-mir-21 hg19:17:57918633-57918656:+        153                      21.791392            N mature,MIMAT0000076 TCGA-CU-A0YN-11A-11R-A10V-13
# 1338 hsa-mir-21 hg19:17:57918633-57918657:+          3                       0.427282            N mature,MIMAT0000076 TCGA-CU-A0YN-11A-11R-A10V-13

hg19_MIMAT0000076_01_grouped_sums <- hg19_MIMAT0000076 %>%
  filter(
    substring(barcode,14,15) == "01"
  ) %>%
  group_by(barcode) %>%
  summarize(sum = sum(reads_per_million_miRNA_mapped))

range(hg19_MIMAT0000076_01_grouped_sums$sum)
# 27564.75 635255.60

table(substr(hg19_MIMAT0000076_01_grouped_sums$barcode, 14,15))
# 01 409


################################################################################
#-- Legacy hg19 - hsa-miR-21-3p is MIMAT0004494
################################################################################
hg19_MIMAT0004494 <- hg19_isoform_datafile[grep("MIMAT0004494", hg19_isoform_datafile$miRNA_region),]
dim(hg19_MIMAT0004494)
# 12971     7
head(hg19_MIMAT0004494)

hg19_MIMAT0004494_01_grouped_sums <- hg19_MIMAT0004494 %>%
  filter(
    substring(barcode,14,15) == "01"
  ) %>%
  group_by(barcode) %>%
  summarize(sum = sum(reads_per_million_miRNA_mapped))
hg19_MIMAT0004494_01_grouped_sums







################################################################################
#-- Main comparison boxplots: hsa-mir-24 family
################################################################################

pdf("miR-24_family.legacy_vs_harmonized.reordered.v3.20190302.pdf", height = 3., width = 4.5)
par(mar = c(3.5,3.,0.5,0.5), mgp = c(3,0.4,0), tck = -0.02)
boxplot(
  log2(hg19_MIMAT0000079_01_grouped_sums$sum), # 24-1-5p
  log2(hg19_MIMAT0004497_01_grouped_sums$sum), # 24-2-5p
  log2(hg19_MIMAT0000080_01_grouped_sums$sum), # 24-3p
  
  log2(hg38_MIMAT0000079_01_grouped_sums$sum), # 24-1-5p
  log2(hg38_MIMAT0004497_01_grouped_sums$sum), # 24-2-5p
  log2(hg38_MIMAT0000080_01_grouped_sums$sum), # 24-3p
  
  log2(hg19_hsa.mir.24.1_stemloop_01_RPMs), # 21-1 stem-loop
  log2(hg19_hsa.mir.24.2_stemloop_01_RPMs), # 21-2 stem-loop
  
  log2(hg38_hsa.mir.24.1_stemloop_01_RPMs), # 21-1 stem-loop
  log2(hg38_hsa.mir.24.2_stemloop_01_RPMs), # 21-2 stem-loop
  
  xlim=c(0.75,10.25),
  ylim=c(-1,16),
  col = c("orange", "skyblue2", "red", "orange", "skyblue2", "red", "grey90", "grey90", "grey", "grey"),
  names = c("24-1-5p", "24-2-5p", "24-3p", "24-1-5p", "24-2-5p", "24-3p", "24-1 SL", "24-2 SL", "24-1 SL", "24-2 SL"),
  cex.axis = 0.8,
  at = 1:10,
  notch = T,
  las = 2,
  varwidth = T,
  cex = 0.6
)
grid(nx = NA, ny = NULL, lty = 1, lwd = 0.5)

lines(c(3.5,3.5), c(-2,17), col = "grey", lwd = 1, lty = 1)
lines(c(6.5,6.5), c(-2,17), col = "grey", lwd = 1, lty = 1)
lines(c(8.5,8.5), c(-2,17), col = "grey", lwd = 1, lty = 1)

mtext("log2(RPM)", side = 2, line = 1.5, cex = 1.2)
#mtext("mir-24 family, BLCA", side = 3, line = 0.2, cex = 1.2)
dev.off()


#-- Wilcoxon pairwise tests
#-- MIMAT0000079: 24-1-5p
wilcox.test(hg19_MIMAT0000079_01_grouped_sums$sum, hg38_MIMAT0000079_01_grouped_sums$sum)$p.value
# 0.4559473

#-- MIMAT0004497: 24-2-5p
wilcox.test(hg19_MIMAT0004497_01_grouped_sums$sum, hg38_MIMAT0004497_01_grouped_sums$sum)$p.value
# 0.5553388

#-- MIMAT0000080 - 24-3p
wilcox.test(hg19_MIMAT0000080_01_grouped_sums$sum, hg38_MIMAT0000080_01_grouped_sums$sum)$p.value
# 0.8537448

#-- mir-24-1 stem-loop
wilcox.test(hg19_hsa.mir.24.1_stemloop_01_RPMs, hg38_hsa.mir.24.1_stemloop_01_RPMs)$p.value
# 1.444217e-136

#-- mir-24-2 stem-loop
wilcox.test(hg19_hsa.mir.24.2_stemloop_01_RPMs, hg38_hsa.mir.24.2_stemloop_01_RPMs)$p.value
# 1.300022e-54
#-- hg38: mir-24-1 and 24-2 stem-loops

wilcox.test(hg38_hsa.mir.24.1_stemloop_01_RPMs, hg38_hsa.mir.24.2_stemloop_01_RPMs)$p.value
# 0.8932015




################################################################################
#-- Main comparison boxplots: hsa-mir-21
################################################################################
pdf("miR-21_family.legacy_vs_harmonized.reordered.grey.20190302.pdf", height = 3., width = 3.5)
par(mar = c(3.5,3.,0.5,0.5), mgp = c(3,0.4,0), tck = -0.02)
rect()
boxplot(
  log2(hg19_MIMAT0000076_01_grouped_sums$sum), # 21-5p
  log2(hg19_MIMAT0004494_01_grouped_sums$sum), # 21-3p
  
  log2(hg38_MIMAT0000076_01_grouped_sums$sum), # 21-5p
  log2(hg38_MIMAT0004494_01_grouped_sums$sum), # 21-3p
  
  log2(hg19_hsa.mir.21_stemloop_RPMs_01), # 21 stem-loop
  log2(hg38_hsa.mir.21_RPMs_stemloop_01), # 21 stem-loop
  
  xlim=c(0.5,6.5),  ylim=c(6,21),
  
  col = c("orange", "red", "orange", "red", "grey90", "grey"),
  names = c("21-5p", "21-3p", "21-5p", "21-3p", "21 SL", "21 SL"),
  cex.axis = 0.8,
  at = 1:6,
  notch = T,
  las = 2,
  varwidth = T,
  cex = 0.6
)
grid(nx = NA, ny = NULL, lty = 1, lwd = 0.5)

lines(c(2.5,2.5), c(-2,22), col = "grey", lwd = 1, lty = 1)
lines(c(4.5,4.5), c(-2,22), col = "grey", lwd = 1, lty = 1)

mtext("log2(RPM)", side = 2, line = 1.5, cex = 1.2)
#mtext("mir-24 family, BLCA", side = 3, line = 0.2, cex = 1.2)
par(bg = 'white')
dev.off()


#-- Wilcoxon pairwise tests
#-- miR-21-5p
wilcox.test(hg19_MIMAT0000076_01_grouped_sums$sum, hg38_MIMAT0000076_01_grouped_sums$sum)$p.value
# 0.8457443

#-- miR-21-3p
wilcox.test(hg19_MIMAT0004494_01_grouped_sums$sum, hg38_MIMAT0004494_01_grouped_sums$sum)$p.value
# 0.8603853

#-- mir-21 stem-loop
wilcox.test(hg19_hsa.mir.21_stemloop_RPMs_01, hg38_hsa.mir.21_RPMs_stemloop_01)$p.value
# 0.8546602

