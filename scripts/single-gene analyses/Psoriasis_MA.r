library(easypackages)
libraries("tidyverse",
          "clusterProfiler",
          "readxl",
          "MetaVolcanoR")
res_SRP055813 <- read_excel("~/Psobesity/SRP055813.txt/Individual_results/res_SRP055813.xlsx")
res_SRP055813 <- as.data.frame(res_SRP055813[!duplicated(res_SRP055813$GeneID), ])
row.names(res_SRP055813) <- res_SRP055813$GeneID
res_SRP055813 <- res_SRP055813[!is.na(res_SRP055813$padj), ]
res_SRP065812 <- read_excel("~/Psobesity/SRP065812.txt/Individual_results/res_SRP065812.xlsx")
res_SRP065812 <- as.data.frame(res_SRP065812[!duplicated(res_SRP065812$GeneID), ])
row.names(res_SRP065812) <- res_SRP065812$GeneID
res_SRP065812 <- res_SRP065812[!is.na(res_SRP065812$padj), ]
res_SRP165679 <- read_excel("~/Psobesity/SRP165679.txt/Individual_results/res_SRP165679.xlsx")
res_SRP165679 <- as.data.frame(res_SRP165679[!duplicated(res_SRP165679$GeneID), ])
row.names(res_SRP165679) <- res_SRP165679$GeneID
res_SRP165679 <- res_SRP165679[!is.na(res_SRP165679$padj), ]
res_SRP145260 <- read_excel("~/Psobesity/SRP145260.txt/Individual_results/res_SRP145260.xlsx")
res_SRP145260 <- as.data.frame(res_SRP145260[!duplicated(res_SRP145260$GeneID), ])
row.names(res_SRP145260) <- res_SRP145260$GeneID
res_SRP145260 <- res_SRP145260[!is.na(res_SRP145260$padj), ]
list_Psoriasis <- list(
  SRP055813 = res_SRP055813,
  SRP065812 = res_SRP065812,
  SRP165679 = res_SRP165679,
  SRP145260 = res_SRP145260
)
common_genes <- Reduce(intersect, lapply(list_Psoriasis, rownames))
for(i in names(list_Psoriasis)){
  list_Psoriasis[[i]] <- as.data.frame(list_Psoriasis[[i]][list_Psoriasis[[i]]$GeneID %in% common_genes, ])
}
setwd("~/Psobesity/Psoriasis")
Psoriasis_meta <- combining_mv(diffexp = list_Psoriasis,
                               pcriteria = "pvalue",
                               genenamecol = "GeneID",
                               metafc = "Mean",
                               foldchangecol = "log2FoldChange",
                               metathr = 0.01,
                               jobname = "MetaVolcano",
                               collaps = TRUE)
res_psoriasis <- Psoriasis_meta@metaresult
res_psoriasis$metap <- p.adjust(res_psoriasis$metap, method = "BH")
genes <- res_SRP055813[, c("GeneID", "Entrez", "Symbol", "Description")]
res_psoriasis <- merge(res_psoriasis, genes, by = "GeneID")
#####GSEA
geneList_psoriasis <- res_psoriasis$Entrez
fold_geneList_psoriasis <- res_psoriasis$metafc
names(fold_geneList_psoriasis) <- as.character(geneList_psoriasis)
fold_geneList_psoriasis <- sort(fold_geneList_psoriasis, decreasing = T)
GSE_psoriasis <- gseKEGG(fold_geneList_psoriasis,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH",
                         organism = 'hsa', 
                         verbose = FALSE
)
#####Save#####
save.image(file = "~/Psobesity/Psoriasis/psoriasis_meta-analysis.RData")
rm(list=ls())
gc()