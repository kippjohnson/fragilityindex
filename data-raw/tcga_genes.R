# data downloaded from http://www.cbioportal.org/study?id=coadread_tcga_pub#summary
# using the GUI option to "download data" on November 7, 2016. From this file, we use
# the table entitles "data_clinical.txt" for survival information. For the CNA patient
# survival data, we use the GUI on the same website to check off genes BCL2L1 and
# POFUT1 and then used the GUI option to download this data, saving the file as
# "CNA_genes_txt".

library(data.table)
library(stringr)

# Read in and prepare Copy Number Alteration TCGA table
in2 <- fread("data-raw/CNAGenes.txt")

in.sub <- in2[Gene %in% c("BCL2L1","POFUT1"),]

bcl.cases <- in2[Gene=="BCL2L1",Cases]
pofut.cases <- in2[Gene=="POFUT1",Cases]

bcl.case0 <- unlist(str_split(bcl.cases,","))
pofut.case0 <- unlist(str_split(pofut.cases,","))

bcl.case1 <- str_sub(bcl.case0, end = -4L)
pofut.case1 <- str_sub(pofut.case0, end = -4L)

# Read in and prepare survival data/clinical TCGA information
in1 <- fread("data-raw/data_clinical.txt")

in1.sub1 <- in1[PATIENT_ID %in% c(bcl.case1, pofut.case1)]
in1.sub2 <- in1[!(PATIENT_ID %in% c(bcl.case1, pofut.case1))]

in1.sub1[,strat:="mutated"]
in1.sub2[,strat:="not_mutated"]

in1.sub1$OS_MONTHS[which(in1.sub1$OS_MONTHS==0)] <- max(in1.sub1$OS_MONTHS)
in1.sub2$OS_MONTHS[which(in1.sub2$OS_MONTHS==0)] <- max(in1.sub2$OS_MONTHS)

in1.sub1$months <- in1.sub1$OS_MONTHS
in1.sub2$months <- in1.sub2$OS_MONTHS

in1.com <- rbind(in1.sub1, in1.sub2)
in1.com$status <- as.factor(in1.com$OS_STATUS)
levels(in1.com$status) <- c(1,0)
in1.com$status <- as.numeric(as.character(in1.com$status))
in1.com$strat <- as.factor(in1.com$strat)

com <- as.data.frame(in1.com)
com <- com[,c("months","status","strat")]

tcga_genes = com

devtools::use_data(tcga_genes, overwrite = TRUE)
