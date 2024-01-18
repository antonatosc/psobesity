#####Libraries#####
library(easypackages)
libraries("tidyverse",
          "readxl",
          "data.table",
          "DESeq2",
          "openxlsx",
          "topconfects",
          "biomaRt")
#####SRP016583#####
countdata <- read.delim("~/Psobesity/SRP055813.txt/final_count", row.names=1, comment.char="#") # nolint
countdata <- countdata[ ,-(1:5)]
metadata <- fread('/media/charis/hdd1/RNA-seq/datasets/SRP055813_metadata')
row.names(metadata) <- metadata$Run
colnames(countdata) <- row.names(metadata)
metadata <- metadata[match(colnames(countdata), metadata$Run), ]
dds_SRP016583 <- DESeqDataSetFromMatrix(countData = countdata,
                                        colData = metadata,
                                        design = ~disease)
keep_SRP016583 <- rowSums(counts(dds_SRP016583)) >= 10
dds_SRP016583 <- dds_SRP016583[keep_SRP016583,]
dds_SRP016583 <- DESeq(dds_SRP016583)
res_SRP016583 <- results(dds_SRP016583)
res_SRP016583 <- as.data.frame(res_SRP016583)
res_SRP016583$GeneID <- rownames(res_SRP016583)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- row.names(res_SRP016583)
G_list <- getBM(filters= "ensembl_gene_id",
  attributes= c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "description"),
  values=genes,
  mart= mart
)
names(G_list) <- c("GeneID", "Entrez", "Symbol", "Description")
res_SRP016583 <- merge(res_SRP016583, G_list, by = "GeneID")
write.xlsx(res_SRP016583, "~/Psobesity/SRP055813.txt/Individual_results/res_SRP055813.xlsx", quote = FALSE)
save.image(file = "~/Psobesity/SRP055813.txt/Individual_results/SRP055813")
rm(list=ls())
gc()
#####SRP065812#####
countdata <- read.delim("~/Psobesity/SRP065812.txt/final_count", row.names=1, comment.char="#")
countdata <- countdata[ ,-(1:5)]
metadata <- fread('/media/charis/hdd1/RNA-seq/datasets/SRP065812_metadata')
row.names(metadata) <- metadata$Run
colnames(countdata) <- row.names(metadata)
metadata <- metadata[match(colnames(countdata), metadata$Run), ]
metadata$subject_status <- sub(" ", "_", metadata$subject_status)
dds_SRP065812 <- DESeqDataSetFromMatrix(countData = countdata,
                                        colData = metadata,
                                        design = ~subject_status)
keep_SRP065812 <- rowSums(counts(dds_SRP065812)) >= 10
dds_SRP065812 <- dds_SRP065812[keep_SRP065812,]
dds_SRP065812 <- DESeq(dds_SRP065812)
res_SRP065812 <- results(dds_SRP065812)
res_SRP065812 <- as.data.frame(res_SRP065812)
res_SRP065812$GeneID <- rownames(res_SRP065812)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- row.names(res_SRP065812)
G_list <- getBM(filters= "ensembl_gene_id",
  attributes= c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "description"),
  values=genes,
  mart= mart
)
names(G_list) <- c("GeneID", "Entrez", "Symbol", "Description")
res_SRP065812 <- merge(res_SRP065812, G_list, by = "GeneID")
write.xlsx(res_SRP065812, "~/Psobesity/SRP065812.txt/Individual_results/res_SRP065812.xlsx", quote = FALSE)
save.image(file = "~/Psobesity/SRP065812.txt/Individual_results/SRP065812")
rm(list=ls())
gc()
#####SRP165679#####
countdata <- read.delim("~/Psobesity/SRP165679.txt/final_count", row.names=1, comment.char="#") # nolint
countdata <- countdata[ ,-(1:5)]
metadata <- fread('/media/charis/hdd1/RNA-seq/datasets/SRP165679_metadata')
row.names(metadata) <- metadata$Run
colnames(countdata) <- row.names(metadata)
metadata <- metadata[match(colnames(countdata), metadata$Run), ]
metadata <- metadata[22:87, ]
countdata <- countdata[, metadata$Run]
dds_SRP165679 <- DESeqDataSetFromMatrix(countData = countdata,
                                        colData = metadata,
                                        design = ~Skin_Type)
keep_SRP165679 <- rowSums(counts(dds_SRP165679)) >= 10
dds_SRP165679 <- dds_SRP165679[keep_SRP165679,]
dds_SRP165679 <- DESeq(dds_SRP165679)
res_SRP165679 <- results(dds_SRP165679)
res_SRP165679 <- as.data.frame(res_SRP165679)
res_SRP165679$GeneID <- rownames(res_SRP165679)
res_SRP165679$UP_CI <- (res_SRP165679$log2FoldChange + res_SRP165679$lfcSE)
res_SRP165679$LO_CI <- (res_SRP165679$log2FoldChange - res_SRP165679$lfcSE)
#####PCA#####
rlog_SRP165679 <- vst(dds_SRP165679, blind = FALSE)
pdf(file = "~/Psobesity/SRP165679.txt/Individual_results/PCA_SRP165679.pdf",
    width = 15,
    height = 10)
plotPCA(rlog_SRP165679, intgroup = "Skin_Type") +
  theme_bw() + 
  geom_point(size = 5) +
  ggtitle(label = "Principal Component Analysis (PCA)"
  ) +
  geom_text(aes(label=name),vjust=2)
dev.off()
#####Save_results#####
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- row.names(res_SRP165679)
G_list <- getBM(filters= "ensembl_gene_id",
                attributes= c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "description"),
                values=genes,
                mart= mart
)
names(G_list) <- c("GeneID", "Entrez", "Symbol", "Description")
res_SRP165679 <- merge(res_SRP165679, G_list, by = "GeneID")
write.xlsx(res_SRP165679, "~/Psobesity/SRP165679.txt/Individual_results/res_SRP165679.xlsx", quote = FALSE)
save.image(file = "~/Psobesity/SRP165679.txt/Individual_results/SRP165679")
rm(list=ls())
gc()
#####SRP145260#####
countdata <- read.delim("~/Psobesity/SRP145260.txt/final_count", row.names=1, comment.char="#") # nolint
countdata <- countdata[ ,-(1:5)]
metadata <- fread('/media/charis/hdd1/RNA-seq/datasets/SRP145260_metadata')
row.names(metadata) <- metadata$Run
colnames(countdata) <- row.names(metadata)
metadata <- metadata[match(colnames(countdata), metadata$Run), ]
dds_SRP145260 <- DESeqDataSetFromMatrix(countData = countdata,
                                        colData = metadata,
                                        design = ~source_name)
keep_SRP145260 <- rowSums(counts(dds_SRP145260)) >= 10
dds_SRP145260 <- dds_SRP145260[keep_SRP145260,]
dds_SRP145260 <- DESeq(dds_SRP145260)
res_SRP145260 <- results(dds_SRP145260)
res_SRP145260 <- as.data.frame(res_SRP145260)
res_SRP145260$GeneID <- rownames(res_SRP145260)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- row.names(res_SRP145260)
G_list <- getBM(filters= "ensembl_gene_id",
                attributes= c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "description"),
                values=genes,
                mart= mart
)
names(G_list) <- c("GeneID", "Entrez", "Symbol", "Description")
res_SRP145260 <- merge(res_SRP145260, G_list, by = "GeneID")
write.xlsx(res_SRP145260, "~/Psobesity/SRP145260.txt/Individual_results/res_SRP145260.xlsx", quote = FALSE)
save.image(file = "~/Psobesity/SRP145260.txt/Individual_results/SRP145260")
rm(list=ls())
gc()