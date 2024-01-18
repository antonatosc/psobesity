library(easypackages)
libraries("tidyverse",
          "clusterProfiler",
          "readxl",
          "data.table",
          "WGCNA",
          "biomaRt",
          "pROC"
)
load('~/Psobesity/Consensus/cWGCNA_bwnet.RData')
consMEs = bwnet$multiMEs
moduleLabels = bwnet$colors
moduleColors = labels2colors(moduleLabels)
bwLabels = matchLabels(bwnet$colors, moduleLabels, pThreshold = 1e-7)
bwColors = labels2colors(bwLabels)
consTree = bwnet$dendrograms[[1]]

#####multi-trait correlation#####
metadata_transform_pso <- metadata_psoriasis
metadata_transform_pso$Status[metadata_transform_pso$Status == "Psoriasis"] <- 2
metadata_transform_pso$Status[metadata_transform_pso$Status == "Healthy"] <- 1
metadata_transform_obs <- metadata_Obesity
metadata_transform_obs$Status[metadata_transform_obs$Status == "Obese"] <- 2
metadata_transform_obs$Status[metadata_transform_obs$Status == "Lean"] <- 1
#calculate correlations for each consensus cluster
Traits = vector(mode="list", length = nSets)
row.names(metadata_transform_pso) <- metadata_transform_pso$Run
Traits[[1]] <- list(data = metadata_transform_pso[ , c(1,4)])
Traits[[2]] <- list(data = metadata_transform_obs[ , c(1,4)])
rownames(Traits[[1]]$data) <- Traits[[1]]$data$Run
rownames(Traits[[2]]$data) <- Traits[[2]]$data$Run
Traits[[1]]$data$Run <- NULL
Traits[[2]]$data$Run <- NULL
moduleTraitCor = list()
moduleTraitPvalue = list()
exprSize = checkSets(multiExpr)
for(i in 1:nSets) {
    moduleTraitCor[[i]] = bicor(consMEs[[i]]$data, as.numeric(Traits[[i]]$data$Status), use = "p")
    moduleTraitPvalue[[i]] = corPvalueFisher(moduleTraitCor[[i]], exprSize$nSamples[i])
    }
psoriasis_correlation <- data.frame(
  "Psoriasis Correlation" = moduleTraitCor[1],
  "Psoriasis P" = moduleTraitPvalue[1],
  "moduleLabel" = rownames(moduleTraitCor[[1]])
)
names(psoriasis_correlation) <- c("Psoriasis Correlation", "Psoriasis P", "moduleLabel")
obesity_correlation <- data.frame(
  "Obesity Correlation" = moduleTraitCor[2],
  "Obesity P" = moduleTraitPvalue[2],
  "moduleLabel" = rownames(moduleTraitCor[[2]])
)
names(obesity_correlation) <- c("Obesity Correlation", "Obesity P", "moduleLabel")
#get conservative consensus correlation patterns
consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))
consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))
#negative cors
negative = moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0
consensusCor[negative] = pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative])
consensusPvalue[negative] = pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative])
#positive cors
positive = moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0
consensusCor[positive] = pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive])
consensusPvalue[positive] = pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive])
rownames(consensusCor) <- rownames(moduleTraitCor[[1]])
consensusCor <- as.data.frame(consensusCor)
consensusCor$moduleLabels <- rownames(consensusCor)
rownames(consensusPvalue) <- rownames(moduleTraitPvalue[[1]])
consensusPvalue <- as.data.frame(consensusPvalue)
consensusPvalue$moduleLabels <- rownames(consensusPvalue)
consensus_patterns <- merge(consensusPvalue, consensusCor, by = "moduleLabels")
names(consensus_patterns) <- c("moduleLabel", "Consensus P", "Consensus Correlation")
#merge all correlations in a dataframe
stats <- Reduce(function(x, y) merge(x, y, all=TRUE), list(psoriasis_correlation,
                                                           obesity_correlation,
                                                           consensus_patterns))  
stats$`Consensus P`[is.na(stats$`Consensus P`)] <- 1
stats$`Consensus Correlation`[is.na(stats$`Consensus Correlation`)] <- 0
melted <- stats[, c(1, 2, 4)]
melted <- melt(melted)
significant <- stats[stats$`Consensus P` < 0.05, ]$moduleLabel
#Fig 3a
ggplot(data = melted,
       aes(x = moduleLabel)) +
  geom_col(data = melted[melted$moduleLabel %in% significant, ],
           aes(x = moduleLabel, y = value, fill = variable),
           stat = "identity",
           position = "dodge",
           width = 1) +
  geom_col(data = melted[!melted$moduleLabel %in% significant, ],
           aes(x = moduleLabel, y = value, fill = variable),
           stat = "identity",
           position = "dodge",
           alpha = 0.2,
           width = 0.5) +
  scale_fill_manual(values = c("red4", "forestgreen")) +
  theme(legend.position = "bottom") +
  coord_flip() +
  bbplot::bbc_style() +
  ylim(c(-1,1))

#####Get modules and enrichment results#####
module_df <- data.frame(
  GeneID = names(bwnet$colors),
  color = labels2colors(bwnet$colors),
  Module = paste0("ME", labels2colors(bwnet$colors)),
  moduleLabels = paste0("ME", moduleLabels)
  )
module_df <- merge(module_df, G_list, by = "GeneID")
enrichment_list <- list()
for(i in unique(module_df$moduleLabels)){
  enrichment_list[[i]] <- enrichKEGG(module_df[module_df$moduleLabels %in% i, ]$Entrez,
                                                organism = 'hsa')
} #Fig 3c

#####eigengene network#####
Status = vector(mode = "list", length = nSets)
for (set in 1:nSets)
{
Status[[set]] = list(data = as.data.frame(Traits[[set]]$data$Status))
names(Status[[set]]$data) = "Status"
}
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
consMEsC = multiSetMEs(multiExpr, universalColors = moduleColors)
MET = consensusOrderMEs(addTraitToMEs(consMEsC, Status))
pdf(file = '~/Psobesity/Consensus/Eigengene.pdf', height = 15, width = 15)
plotEigengeneNetworks(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
  zlimPreservation = c(0.5, 1), xLabelsAngle = 90) #Fig 3d
dev.off()

#####Network export#####
network_export <- function(module, threshold){
  cons_module <- module_df[module_df$moduleLabels == module, ]
  cons_module <- cons_module[!duplicated(cons_module$GeneID), ]
  exprs_pso <- norm_psoriasis[, cons_module$GeneID]
  exprs_obs <- norm_Obesity[, cons_module$GeneID]
  multiExpr = vector(mode = "list", length = nSets)
  multiExpr[[1]] = list(data = as.data.frame(exprs_pso))
  names(multiExpr[[1]]$data) <- colnames(exprs_pso)
  rownames(multiExpr[[1]]$data) <- rownames(exprs_pso)
  multiExpr[[2]] = list(data = as.data.frame(exprs_obs))
  names(multiExpr[[2]]$data) <- colnames(exprs_obs)
  rownames(multiExpr[[2]]$data) <- rownames(exprs_obs)
  cor <- WGCNA::bicor
  #we reconstruct TOM matrix according to the chosen module
  adjacencies = array(0, dim = c(nSets, length(colnames(multiExpr[[1]]$data)), length(colnames(multiExpr[[1]]$data))))
  TOM = array(0, dim = c(nSets, length(colnames(multiExpr[[1]]$data)), length(colnames(multiExpr[[1]]$data))))
  for (i in 1:nSets) {
    #adjacencies in each set
    adjacencies[i, , ] = abs(WGCNA::bicor(multiExpr[[i]]$data, use = "p"))^12
    #TOM matrix in each set
    TOM[i, , ] = TOMsimilarity(adjacencies[i, , ],
                              TOMType = "signed")
  }
  consensusTOM = pmin(TOM[1, , ], TOM[2, , ])
  probes = names(multiExpr[[1]]$data)
  dimnames(consensusTOM) = list(probes, probes)
  names_probes <- module_df[module_df$GeneID %in% probes, ]
  names_probes <- names_probes[!duplicated(names_probes$GeneID), ]
  cyt = exportNetworkToCytoscape(consensusTOM,
                                 weighted = TRUE,
                                 threshold = threshold,
                                 nodeNames = probes,
                                 altNodeNames = names_probes$Symbol)
  return(cyt)
} #Visualize in cytoscape for Fig 3e

#####ROC#####
roc_analysis <- function(disease){
  ifelse(disease == "psoriasis",
         number <- 1,
         number <- 2)
  ifelse(disease == "psoriasis",
         metadata <- metadata_psoriasis,
         metadata <- metadata_Obesity)
  dataset <- as.data.frame(consMEs[[number]])
  dataset$Run <- rownames(dataset)
  dataset <- merge(metadata, dataset)
  dataset_roc <- list()
  results <- data.frame()
  for(i in paste0("data.ME", 0:50)){
    dataset_roc[[i]] <- roc(dataset$Status,
                            as.numeric(dataset[[i]]),
                            ci = TRUE,
                            smooth = TRUE,
                            n.bootstrap = 2000)
    results <- rbind(results, data.frame(
      "C-statistic" = dataset_roc[[i]]$ci[2],
      "95% CI" = paste0(dataset_roc[[i]]$ci[1], "-", dataset_roc[[i]]$ci[3]),
      "Module" = i))
  }
    return(results)
  }

#####Dendrogram#####
plotDendroAndColors(
  consTree,
  moduleColors,
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )
  module_df <- data.frame(
  gene_id = names(bwnet$colors),
  colors = labels2colors(bwnet$colors)
)