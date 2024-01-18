library(easypackages)
libraries("tidyverse",
          "clusterProfiler",
          "readxl",
          "MetaVolcanoR")
res_SRP132990 <- read_excel("~/Psobesity/SRP132990.txt/Individual_results/res_SRP132990.xlsx")
res_SRP132990 <- as.data.frame(res_SRP132990[!duplicated(res_SRP132990$GeneID), ])
row.names(res_SRP132990) <- res_SRP132990$GeneID
res_SRP132990 <- res_SRP132990[!is.na(res_SRP132990$padj), ]
res_SRP268322 <- read_excel("~/Psobesity/SRP268322.txt/Individual_results/res_SRP268322.xlsx")
res_SRP268322 <- as.data.frame(res_SRP268322[!duplicated(res_SRP268322$GeneID), ])
row.names(res_SRP268322) <- res_SRP268322$GeneID
res_SRP268322 <- res_SRP268322[!is.na(res_SRP268322$padj), ]
res_SRP278883 <- read_excel("~/Psobesity/SRP278883.txt/Individual_results/res_SRP278883.xlsx")
res_SRP278883 <- as.data.frame(res_SRP278883[!duplicated(res_SRP278883$GeneID), ])
row.names(res_SRP278883) <- res_SRP278883$GeneID
res_SRP278883 <- res_SRP278883[!is.na(res_SRP278883$padj), ]
res_SRP295864 <- read_excel("~/Psobesity/SRP295864.txt/Individual_results/res_SRP295864.xlsx")
res_SRP295864 <- as.data.frame(res_SRP295864[!duplicated(res_SRP295864$GeneID), ])
row.names(res_SRP295864) <- res_SRP295864$GeneID
res_SRP295864 <- res_SRP295864[!is.na(res_SRP295864$padj), ]
res_SRP304398 <- read_excel("~/Psobesity/SRP304398.txt/Individual_results/res_SRP304398.xlsx")
res_SRP304398 <- as.data.frame(res_SRP304398[!duplicated(res_SRP304398$GeneID), ])
row.names(res_SRP304398) <- res_SRP304398$GeneID
res_SRP304398 <- res_SRP304398[!is.na(res_SRP304398$padj), ]
list_Obesity <- list(
    SRP132990 = res_SRP132990,
    SRP268322 = res_SRP268322,
    SRP278883 = res_SRP278883,
    SRP295864 = res_SRP295864,
    SRP304398 = res_SRP304398
)
common_genes <- Reduce(intersect, lapply(list_Obesity, rownames))
for(i in names(list_Obesity)){
  list_Obesity[[i]] <- as.data.frame(list_Obesity[[i]][list_Obesity[[i]]$GeneID %in% common_genes, ])
}
setwd("~/Psobesity/Obesity")
Obesity_meta <- combining_mv(diffexp = list_Obesity,
                             pcriteria = "pvalue",
                             genenamecol = "GeneID",
                             metafc = "Mean",
                             foldchangecol = "log2FoldChange",
                             metathr = 0.01,
                             jobname = "MetaVolcano",
                             collaps = TRUE)
res_Obesity <- Obesity_meta@metaresult
res_Obesity$metap <- p.adjust(res_Obesity$metap, method = "BH")
genes <- res_SRP132990[, c("GeneID", "Entrez", "Symbol", "Description")]
res_Obesity <- merge(res_Obesity, genes, by = "GeneID")
#####GSEA
geneList_Obesity <- res_Obesity$Entrez
fold_geneList_Obesity <- res_Obesity$randomSummary
names(fold_geneList_Obesity) <- as.character(geneList_Obesity)
fold_geneList_Obesity <- sort(fold_geneList_Obesity, decreasing = TRUE)
GSE_Obesity <- gseKEGG(fold_geneList_Obesity, 
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH",
                         organism = 'hsa', 
                         verbose = FALSE
                )
#####Save#####
save.image(file = "~/Psobesity/Obesity/Obesity")
rm(list=ls())
gc()