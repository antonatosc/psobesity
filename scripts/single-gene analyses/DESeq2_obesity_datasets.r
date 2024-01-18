#####Libraries#####
library(easypackages)
libraries("tidyverse",
          "readxl",
          "data.table",
          "DESeq2",
          "openxlsx",
          "topconfects",
          "biomaRt")
#####SRP132990#####
countdata <- read.delim("~/Psobesity/SRP132990.txt/final_count", row.names=1, comment.char="#")
countdata <- countdata[ ,-(1:5)]
metadata <- fread('~/Psobesity/datasets/SRP132990_metadata')
row.names(metadata) <- metadata$Run
colnames(countdata) <- row.names(metadata)
metadata <- metadata[match(colnames(countdata), metadata$Run), ]
dds_SRP132990 <- DESeqDataSetFromMatrix(countData = countdata,
                                        colData = metadata,
                                        design = ~disease_state)
dds_SRP132990 <- collapseReplicates(dds_SRP132990,
                             dds_SRP132990$BioSample,
                             dds_SRP132990$Run)
keep_SRP132990 <- rowSums(counts(dds_SRP132990)) >= 10
dds_SRP132990 <- dds_SRP132990[keep_SRP132990,]
dds_SRP132990 <- DESeq(dds_SRP132990)
res_SRP132990 <- results(dds_SRP132990)
res_SRP132990 <- as.data.frame(res_SRP132990)
res_SRP132990$GeneID <- rownames(res_SRP132990)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- row.names(res_SRP132990)
G_list <- getBM(filters= "ensembl_gene_id",
                attributes= c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "description"),
                values=genes,
                mart= mart
)
names(G_list) <- c("GeneID", "Entrez", "Symbol", "Description")
res_SRP132990 <- merge(res_SRP132990, G_list, by = "GeneID")
write.xlsx(res_SRP132990, "~/Psobesity/SRP132990.txt/Individual_results/res_SRP132990.xlsx", quote = FALSE)
save.image(file = "~/Psobesity/SRP132990.txt/Individual_results/SRP132990.RData")
rm(list=ls())
gc()
#####SRP304398#####
countdata <- read.delim("~/Psobesity/SRP304398.txt/final_count", row.names=1, comment.char="#")
countdata <- countdata[ ,-(1:5)]
metadata <- fread('~/Psobesity/datasets/SRP304398_metadata')
row.names(metadata) <- metadata$Run
colnames(countdata) <- row.names(metadata)
metadata <- metadata[match(colnames(countdata), metadata$Run), ]
metadata$disease <- c("lean", "lean", "lean", "lean", "lean", "obese", "obese", "obese", "obese", "obese", "obese", "obese", "obese", "obese", "obese", "obese", "obese", "obese", "obese", "obese")
metadata <- metadata[1:15, ]
countdata <- countdata[, metadata$Run]
dds_SRP304398 <- DESeqDataSetFromMatrix(countData = countdata,
                                        colData = metadata,
                                        design = ~disease)
keep_SRP304398 <- rowSums(counts(dds_SRP304398)) >= 10
dds_SRP304398 <- dds_SRP304398[keep_SRP304398,]
dds_SRP304398 <- DESeq(dds_SRP304398)
res_SRP304398 <- results(dds_SRP304398)
res_SRP304398 <- as.data.frame(res_SRP304398)
res_SRP304398$GeneID <- rownames(res_SRP304398)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- row.names(res_SRP304398)
G_list <- getBM(filters= "ensembl_gene_id",
                attributes= c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "description"),
                values=genes,
                mart= mart
)
names(G_list) <- c("GeneID", "Entrez", "Symbol", "Description")
res_SRP304398 <- merge(res_SRP304398, G_list, by = "GeneID")
write.xlsx(res_SRP304398, "~/Psobesity/SRP304398.txt/Individual_results/res_SRP304398.xlsx", quote = FALSE)
save.image(file = "~/Psobesity/SRP304398.txt/Individual_results/SRP304398.RData")
rm(list=ls())
gc()
#####SRP268322#####
countdata <- read.delim("~/Psobesity/SRP268322.txt/final_count", row.names=1, comment.char="#") # nolint
countdata <- countdata[ ,-(1:5)]
metadata <- fread('~/Psobesity/datasets/SRP268322_metadata')
row.names(metadata) <- metadata$Run
colnames(countdata) <- row.names(metadata)
metadata <- metadata[match(colnames(countdata), metadata$Run), ]
asd <- metadata
dds_SRP268322 <- DESeqDataSetFromMatrix(countData = countdata,
                                        colData = metadata,
                                        design = ~subject_status)
dds_SRP268322 <- collapseReplicates(dds_SRP268322,
                             dds_SRP268322$BioSample,
                             dds_SRP268322$Run)
keep_SRP268322 <- rowSums(counts(dds_SRP268322)) >= 10
dds_SRP268322 <- dds_SRP268322[keep_SRP268322,]
dds_SRP268322 <- DESeq(dds_SRP268322)
res_SRP268322 <- results(dds_SRP268322)
res_SRP268322 <- as.data.frame(res_SRP268322)
res_SRP268322$GeneID <- rownames(res_SRP268322)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- row.names(res_SRP268322)
G_list <- getBM(filters= "ensembl_gene_id",
                attributes= c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "description"),
                values=genes,
                mart= mart
)
names(G_list) <- c("GeneID", "Entrez", "Symbol", "Description")
res_SRP268322 <- merge(res_SRP268322, G_list, by = "GeneID")
write.xlsx(res_SRP268322, "~/Psobesity/SRP268322.txt/Individual_results/res_SRP268322.xlsx", quote = FALSE)
save.image(file = "~/Psobesity/SRP268322.txt/Individual_results/SRP268322.RData")
rm(list=ls())
gc()
#####SRP295864#####
countdata <- read.delim("~/Psobesity/SRP295864.txt/final_count", row.names=1, comment.char="#") # nolint
countdata <- countdata[ ,-(1:5)]
metadata <- fread('~/Psobesity/datasets/SRP295864_metadata')
row.names(metadata) <- metadata$Run
colnames(countdata) <- row.names(metadata)
metadata <- metadata[match(colnames(countdata), metadata$Run), ]
dds_SRP295864 <- DESeqDataSetFromMatrix(countData = countdata,
                                        colData = metadata,
                                        design = ~bmi_group)
keep_SRP295864 <- rowSums(counts(dds_SRP295864)) >= 10
dds_SRP295864 <- dds_SRP295864[keep_SRP295864,]
dds_SRP295864 <- DESeq(dds_SRP295864)
res_SRP295864 <- results(dds_SRP295864)
res_SRP295864 <- as.data.frame(res_SRP295864)
res_SRP295864$GeneID <- rownames(res_SRP295864)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- row.names(res_SRP295864)
G_list <- getBM(filters= "ensembl_gene_id",
                attributes= c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "description"),
                values=genes,
                mart= mart
)
names(G_list) <- c("GeneID", "Entrez", "Symbol", "Description")
res_SRP295864 <- merge(res_SRP295864, G_list, by = "GeneID")
write.xlsx(res_SRP295864, "~/Psobesity/SRP295864.txt/Individual_results/res_SRP295864.xlsx", quote = FALSE)
save.image(file = "~/Psobesity/SRP295864.txt/Individual_results/SRP295864.RData")
rm(list=ls())
gc()
#####SRP278883#####
countdata <- read.delim("~/Psobesity/SRP278883.txt/final_count", row.names=1, comment.char="#") # nolint
countdata <- countdata[ ,-(1:5)]
metadata <- fread('~/Psobesity/datasets/SRP278883_metadata')
row.names(metadata) <- metadata$Run
colnames(countdata) <- row.names(metadata)
metadata <- metadata[match(colnames(countdata), metadata$Run), ]
metadata <- metadata[!metadata$Condition == "MUO", ]
countdata <- countdata[, colnames(countdata) == metadata$Run]
metadata <- metadata[match(colnames(countdata), metadata$Run), ]
dds_SRP278883 <- DESeqDataSetFromMatrix(countData = countdata,
                                        colData = metadata,
                                        design = ~Condition)
keep_SRP278883 <- rowSums(counts(dds_SRP278883)) >= 10
dds_SRP278883 <- dds_SRP278883[keep_SRP278883,]
dds_SRP278883 <- DESeq(dds_SRP278883)
res_SRP278883 <- results(dds_SRP278883)
res_SRP278883 <- as.data.frame(res_SRP278883)
res_SRP278883$GeneID <- rownames(res_SRP278883)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- row.names(res_SRP278883)
G_list <- getBM(filters= "ensembl_gene_id",
                attributes= c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "description"),
                values=genes,
                mart= mart
)
res_SRP278883$GeneID <- genes
names(G_list) <- c("GeneID", "Entrez", "Symbol", "Description")
res_SRP278883 <- merge(res_SRP278883, G_list, by = "GeneID")
write.xlsx(res_SRP278883, "~/Psobesity/SRP278883.txt/Individual_results/res_SRP278883.xlsx", quote = FALSE)
save.image(file = "~/Psobesity/SRP278883.txt/Individual_results/SRP278883.RData")
rm(list=ls())
gc()