#####Consensus#####
library(easypackages)
libraries("tidyverse",
          "readxl",
          "data.table",
          "DESeq2",
          "openxlsx",
          "org.Hs.eg.db",
          "topconfects",
          "apeglm",
          "pheatmap",
          "biomaRt",
          "viridis",
          "gplots",
          "sva",
          "gg3D",
          "genefilter",
          "WGCNA",
          "biomaRt",
          "igraph"
)
#####Obesity
#SRP278883
countdata_SRP278883 <- read.delim("~/Psobesity/SRP278883.txt/final_count", row.names=1, comment.char="#")
countdata_SRP278883 <- countdata_SRP278883[ ,-(1:5)]
metadata_SRP278883 <- fread('~/Psobesity/datasets/SRP278883_metadata')
row.names(metadata_SRP278883) <- metadata_SRP278883$Run
colnames(countdata_SRP278883) <- row.names(metadata_SRP278883)
metadata_SRP278883 <- metadata_SRP278883[match(colnames(countdata_SRP278883), metadata_SRP278883$Run), ]
metadata_SRP278883 <- metadata_SRP278883[!metadata_SRP278883$Condition == "MUO", ]
countdata_SRP278883 <- countdata_SRP278883[, colnames(countdata_SRP278883) == metadata_SRP278883$Run]
metadata_SRP278883 <- metadata_SRP278883[match(colnames(countdata_SRP278883), metadata_SRP278883$Run), ]
names(metadata_SRP278883)[names(metadata_SRP278883) == "Condition"] <- "Status"
metadata_3 <- metadata_SRP278883[, c(1, 3, 4, 9, 16, 17)]
metadata_3$Study <- "SRP278883"
metadata_3$Status <- replace(metadata_3$Status, metadata_3$Status == "MHL", "Lean")
metadata_3$Status <- replace(metadata_3$Status, metadata_3$Status == "MHO", "Obese")
#SRP304398
countdata_SRP304398 <- read.delim("~/Psobesity/SRP304398.txt/final_count", row.names=1, comment.char="#")
countdata_SRP304398 <- countdata_SRP304398[ ,-(1:5)]
metadata_SRP304398 <- fread('~/Psobesity/datasets/SRP304398_metadata')
row.names(metadata_SRP304398) <- metadata_SRP304398$Run
colnames(countdata_SRP304398) <- row.names(metadata_SRP304398)
metadata_SRP304398 <- metadata_SRP304398[match(colnames(countdata_SRP304398), metadata_SRP304398$Run), ]
metadata_SRP304398$disease <- c("lean", "lean", "lean", "lean", "lean", "obese", "obese", "obese", "obese", "obese", "obese", "obese", "obese", "obese", "obese", "obese", "obese", "obese", "obese", "obese")
names(metadata_SRP304398)[names(metadata_SRP304398) == "disease"] <- "Status"
metadata_5 <- metadata_SRP304398[, c(1, 4, 5, 28, 18, 19)]
metadata_5$Study <- "SRP304398"
metadata_5$Status <- replace(metadata_5$Status, metadata_5$Status == "lean", "Lean")
metadata_5$Status <- replace(metadata_5$Status, metadata_5$Status == "obese", "Obese")
metadata_5 <- metadata_5[1:15, ]
countdata_SRP304398 <- countdata_SRP304398[, metadata_5$Run]
#combined
metadata_Obesity <- rbind(metadata_3,
                          metadata_5)
metadata_Obesity<- metadata_Obesity %>% mutate_all( ~ str_replace_all(., " ", "_"))
countdata_obesity <- cbind(countdata_SRP304398,
                           countdata_SRP278883)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- row.names(countdata_obesity)
G_list <- getBM(filters= "ensembl_gene_id",
  attributes= c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "description"),
  values=genes,
  mart= mart
)
names(G_list) <- c("GeneID", "Entrez", "Symbol", "Description") #get gene names
countdata_obesity <- countdata_obesity[, colnames(countdata_obesity) %in% metadata_Obesity$Run]
countdata_obesity$variance = apply(countdata_obesity, 1, var)
countdata_obesity = countdata_obesity[countdata_obesity$variance >= quantile(countdata_obesity$variance, c(.25)), ] #variance filtering
countdata_obesity$variance <- NULL
row.names(metadata_Obesity) <- metadata_Obesity$Run
reorder <- match(colnames(countdata_obesity), rownames(metadata_Obesity))
metadata_Obesity <- metadata_Obesity[reorder, ]
dds_obesity <- DESeqDataSetFromMatrix(countData = countdata_obesity,
                                      colData = metadata_Obesity,
                                      design = ~Status)
vst_obesity <- vst(dds_obesity)
adjusted_count_obesity <- ComBat_seq(as.matrix(countdata_obesity),
                                     batch = metadata_Obesity$Study,
                                     group = metadata_Obesity$Status) #Adjust for Study batch
adjusted_dds_obesity <- DESeqDataSetFromMatrix(countData = adjusted_count_obesity,
                                      colData = metadata_Obesity,
                                      design = ~Status)
adjusted_vst_obesity <- vst(adjusted_dds_obesity)
vst_Obesity <- vst(adjusted_dds_obesity)
norm_Obesity <- as.data.frame(assay(vst_Obesity)) %>%
    t()
pdf(file = "~/Psobesity/Consensus/PCA_obesity.pdf",
    width = 15,
    height = 8)
pca_data <- plotPCA(vst_obesity,
                    intgroup = c("Status", "Study"),
                    ntop = 1000,
                    returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
ggplot(pca_data, aes(PC1, PC2, color=Status, shape=Study)) +
  theme_classic() +
  scale_color_manual(values = c("green", "red")) +
  geom_point(size = 2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Principal Component Analysis before batch effect removal")
adjusted_pca_data <- plotPCA(adjusted_vst_obesity,
                      intgroup = c("Status", "Study"),
                      ntop = 500,
                      returnData=TRUE)
adjusted_percentVar <- round(100 * attr(adjusted_pca_data, "percentVar"))
ggplot(adjusted_pca_data, aes(PC1, PC2, color=Status, shape=Study)) +
  theme_classic() +
  scale_color_manual(values = c("green", "red")) +
  geom_point(size = 2) +
  xlab(paste0("PC1: ",adjusted_percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",adjusted_percentVar[2],"% variance")) +
  ggtitle("Principal Component Analysis after batch effect removal")
dev.off()
#####Psoriasis
#SRP065812
countdata_SRP065812 <- read.delim("~/Psobesity/SRP065812.txt/final_count", row.names=1, comment.char="#")
countdata_SRP065812 <- countdata_SRP065812[ ,-(1:5)]
metadata_SRP065812 <- fread('~/Psobesity/datasets/SRP065812_metadata')
row.names(metadata_SRP065812) <- metadata_SRP065812$Run
colnames(countdata_SRP065812) <- row.names(metadata_SRP065812)
metadata_SRP065812 <- metadata_SRP065812[match(colnames(countdata_SRP065812), metadata_SRP065812$Run), ]
names(metadata_SRP065812)[names(metadata_SRP065812) == "source_name"] <- "Status"
metadata_2 <- metadata_SRP065812[, c(1, 3, 4, 23, 15, 16)]
metadata_2$Study <- "SRP065812"
metadata_2$Status <- replace(metadata_2$Status, metadata_2$Status == "Pre-Treated Psoriasis", "Psoriasis")
metadata_2$Status <- replace(metadata_2$Status, metadata_2$Status == "Normal Skin", "Healthy")
#SRP165679
countdata_SRP165679 <- read.delim("~/Psobesity/SRP165679.txt/final_count", row.names=1, comment.char="#")
countdata_SRP165679 <- countdata_SRP165679[ ,-(1:5)]
metadata_SRP165679 <- fread('~/Psobesity/datasets/SRP165679_metadata')
row.names(metadata_SRP165679) <- metadata_SRP165679$Run
colnames(countdata_SRP165679) <- row.names(metadata_SRP165679)
metadata_SRP165679 <- metadata_SRP165679[match(colnames(countdata_SRP165679), metadata_SRP165679$Run), ]
names(metadata_SRP165679)[names(metadata_SRP165679) == "Skin_Type"] <- "Status"
metadata_3 <- metadata_SRP165679[, c(1, 3, 4, 23, 15, 16)]
metadata_3$Study <- "SRP165679"
metadata_3$Status <- replace(metadata_3$Status, metadata_3$Status == "lesional", "Psoriasis")
metadata_3$Status <- replace(metadata_3$Status, metadata_3$Status == "healthy", "Healthy")
metadata_3 <- metadata_3[22:87, ]
countdata_SRP165679 <- countdata_SRP165679[, metadata_3$Run]
#combined
metadata_psoriasis <- rbind(metadata_2,
                            metadata_3)
metadata_psoriasis<- metadata_psoriasis %>% mutate_all( ~ str_replace_all(., " ", "_"))
countdata_psoriasis <- cbind(countdata_SRP065812,
                             countdata_SRP165679)
countdata_psoriasis <- countdata_psoriasis[, colnames(countdata_psoriasis) %in% metadata_psoriasis$Run]
countdata_psoriasis$variance = apply(countdata_psoriasis, 1, var)
countdata_psoriasis = countdata_psoriasis[countdata_psoriasis$variance >= quantile(countdata_psoriasis$variance, c(.25)), ] #variance filtering 
countdata_psoriasis$variance <- NULL
adjusted_count_psoriasis <- ComBat_seq(as.matrix(countdata_psoriasis),
                                       batch = metadata_psoriasis$Study,
                                       group = metadata_psoriasis$Status) #Adjust for Study batch
dds_psoriasis <- DESeqDataSetFromMatrix(countData = countdata_psoriasis,
                                        colData = metadata_psoriasis,
                                        design = ~Status)
vst_psoriasis <- vst(dds_psoriasis)
adjusted_dds_psoriasis <- DESeqDataSetFromMatrix(countData = adjusted_count_psoriasis,
                                      colData = metadata_psoriasis,
                                      design = ~Status)
adjusted_vst_psoriasis <- vst(adjusted_dds_psoriasis)
pdf(file = "~/Psobesity/Consensus/PCA_psoriasis.pdf",
    width = 15,
    height = 8)
pca_data <- plotPCA(vst_psoriasis,
                    intgroup = c("Status", "Study"),
                    returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
ggplot(pca_data, aes(PC1, PC2, color=Status, shape=Study)) +
  theme_classic() +
  scale_color_manual(values = c("green", "red")) +
  geom_point(size = 2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Principal Component Analysis before batch effect removal")
adjusted_pca_data <- plotPCA(adjusted_vst_psoriasis,
                            intgroup = c("Status", "Study"),
                            returnData=TRUE)
adjusted_percentVar <- round(100 * attr(adjusted_pca_data, "percentVar"))
ggplot(adjusted_pca_data, aes(PC1, PC2, color=Status, shape=Study)) +
  theme_classic() +
  scale_color_manual(values = c("green", "red")) +
  geom_point(size = 2) +
  xlab(paste0("PC1: ",adjusted_percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",adjusted_percentVar[2],"% variance")) +
  ggtitle("Principal Component Analysis after batch effect removal")
dev.off()
vst_psoriasis <- vst(adjusted_dds_psoriasis)
norm_psoriasis <- as.data.frame(assay(vst_psoriasis)) %>%
    t()
enableWGCNAThreads()
#####consensus WGCNA
multiExpr <- list(
    Psoriasis = norm_psoriasis,
    Obesity = norm_Obesity
)
#estimate soft thresholding power for cWGCNA
powers = c(seq(4,10,by=1), seq(12,20, by=2))
nSets = 2
powerTables = vector(mode = "list", length = nSets)
for (set in 1:nSets)
powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]], 
    powerVector=powers,
    corFnc = WGCNA::bicor,
    networkType = "signed",
    verbose = 2)[[2]]
    )
collectGarbage()
colors = c("black", "red")
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
"Max connectivity")
ylim = matrix(NA, nrow = 2, ncol = 4)
setLabels = c("Psoriasis", "Obesity")
for (set in 1:nSets)
{
for (col in 1:length(plotCols))
{
ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
}
}
sizeGrWindow(8, 6)
pdf('~/Psobesity/Consensus/All_Power.pdf', 
height = 10, width = 10)
par(mfcol = c(2,2))
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
if (set==1)
{
plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
main = colNames[col])
addGrid()
}
if (col==1)
{
text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
labels=powers,cex=cex1,col=colors[set])
} else
text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
labels=powers,cex=cex1,col=colors[set])
if (col==1)
{
legend("bottomright", legend = setLabels, col = colors, pch = 20) 
} else
legend("topright", legend = setLabels, col = colors, pch = 20) 
}
dev.off()
setwd('~/Psobesity/Consensus')
multiExpr = vector(mode = "list", length = nSets)
genes <- intersect(colnames(norm_Obesity), colnames(norm_psoriasis))
norm_Obesity <- norm_Obesity[, genes]
norm_psoriasis <- norm_psoriasis[, genes]
multiExpr[[1]] = list(data = as.data.frame(norm_psoriasis))
names(multiExpr[[1]]$data) <- colnames(norm_psoriasis)
rownames(multiExpr[[1]]$data) <- rownames(norm_psoriasis)
multiExpr[[2]] = list(data = as.data.frame(norm_Obesity))
names(multiExpr[[2]]$data) <- colnames(norm_Obesity)
rownames(multiExpr[[2]]$data) <- rownames(norm_Obesity)
cor <- WGCNA::bicor
gc()
bwnet <- blockwiseConsensusModules(multiExpr, 
    maxBlockSize = 50000,
    corType = "bicor",
    networkType = "signed",
    TOMType = "signed",
    power = 12,
    minModuleSize = 20,
    deepSplit = 4,
    pamRespectsDendro = FALSE,
    mergeCutHeight = 0.2,
    numericLabels = TRUE,
    minKMEtoStay = 0,
    saveTOMs = TRUE,
    verbose = 10
    )
rm(list=setdiff(ls(), c("bwnet",
                        "metadata_psoriasis",
                        "metadata_Obesity",
                        "nSets",
                        "multiExpr",
                        "norm_Obesity",
                        "norm_psoriasis",
                        "setLabels",
                        "G_list")))
setwd('~/Psobesity/WGCNA')
save.image('cWGCNA_bwnet.RData')