library(easypackages)
libraries("tidyverse",
          "clusterProfiler",
          "readxl",
          "MetaVolcanoR",
          "ggplot2",
          "cowplot")
setwd("~/Psobesity/Individual_results")
psoriasis_datasets = c("SRP145260", "SRP065812", "SRP055813", "SRP165679")
obesity_datasets = c("SRP278883", "SRP295864", "SRP132990", "SRP268322", "SRP304398")
psoriasis=list()
for (i in psoriasis_datasets) {
  print(paste0("Removing ", i))
  datasets = psoriasis_datasets[psoriasis_datasets != i]
  psoriasis[[i]] = loo("psoriasis", datasets, "obesity", i)
}
obesity = list()
for (i in obesity_datasets) {
  print(paste0("Removing ", i))
  datasets = obesity_datasets[obesity_datasets != i]
  obesity[[i]] = loo("obesity", datasets, "psoriasis", i)
}
gene_overlap_test = function(gene_list1, gene_list2, background, verbose = T){
  A = length(gene_list1)
  B = length(gene_list2)
  C = background
  AB = length(intersect(gene_list1, gene_list2))
  pval = phyper(AB-1, B, C-B, A, lower.tail = FALSE)
  return(pval)
}
loo = function(disease1, datasets, disease2, removed){
  ls = list()
  for(i in datasets){
    print(paste0("Reading:", i))
    ls[[i]] = read_excel(paste0("res_", i, ".xlsx"))
    ls[[i]] = as.data.frame(ls[[i]][!duplicated(ls[[i]][["GeneID"]]), ])
    row.names(ls[[i]]) <- ls[[i]][["GeneID"]]
    ls[[i]] <- ls[[i]][!is.na(ls[[i]][["padj"]]), ]
    common_genes <- Reduce(intersect, lapply(ls, rownames))
    ls[[i]] = as.data.frame(ls[[i]][ls[[i]][["GeneID"]] %in% common_genes, ])
  }
  meta <- combining_mv(diffexp = ls,
                       pcriteria = "pvalue",
                       genenamecol = "GeneID",
                       metafc = "Mean",
                       foldchangecol = "log2FoldChange",
                       metathr = 0.01,
                       jobname = "MetaVolcano",
                       collaps = TRUE)
  res_meta <- meta@metaresult
  res_meta$metap <- p.adjust(res_meta$metap, method = "BH")
  genes <- ls[[1]][, c("GeneID", "Entrez", "Symbol", "Description")]
  res_meta <- merge(res_meta, genes, by = "GeneID")
  #####GSEA
  geneList <- res_meta$Entrez
  fold_geneList <- res_meta$metafc
  names(fold_geneList) <- as.character(geneList)
  fold_geneList <- sort(fold_geneList, decreasing = T)
  GSE1 <- gseKEGG(fold_geneList,
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  organism = 'hsa', 
                  verbose = FALSE
  )
  options(ggrepel.max.overlaps = Inf)
  pathways1 <- GSE1@result
  res_meta2 = read_excel(paste0("res_", disease2, ".xlsx"))
  geneList2 <- res_meta2$Entrez
  fold_geneList2 <- res_meta2$metafc
  names(fold_geneList2) <- as.character(geneList2)
  fold_geneList2 <- sort(fold_geneList2, decreasing = T)
  GSE2 <- gseKEGG(fold_geneList2,
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  organism = 'hsa', 
                  verbose = FALSE
  ) 
  pathways2 <- GSE2@result
  shared1 <- pathways1[pathways1$ID %in% pathways2$ID, ]
  shared2 <- pathways2[pathways2$ID %in% pathways1$ID, ]
  shared <- merge(shared1, shared2, by = "Description")
  p1 = ggplot(shared,
              aes(x = NES.x,
                  y = NES.y,
                  color = ifelse(NES.x>0 & NES.y>0, "blue", 
                                 ifelse(NES.x < 0 & NES.y < 0, "purple",
                                        ifelse(NES.x>0 & NES.y<0, "green", "red"))))) +
    geom_point(size = 4,
               alpha = 0.5) +
    stat_ellipse() +
    geom_point(size = 4,
               alpha = 0.5) +
    bbplot::bbc_style() +
    theme(legend.position = "none",
          axis.title = element_text(size = 18)) +
    xlab(paste0("NES in ", disease1)) +
    ylab(paste0("NES in ", disease2)) +
    theme(text = element_text(size = 15)) +
    ggrepel::geom_label_repel(label = ifelse(shared$Description %in% "Th17 cell differentiation", "Th17 cell differentiation", ""),
                             size = 7,
                             color = "black")
  ####Comparisons
  regulated_1 <- res_meta[(res_meta$metap < 0.01) & !between(res_meta$metafc, -0.58, 0.58), ]
  up_1 = regulated_1[regulated_1$metafc > 0.58, ]
  down_1 = regulated_1[regulated_1$metafc < -0.58, ]
  regulated_2 <- res_meta2[(res_meta2$metap < 0.01) & !between(res_meta2$metafc, -0.58, 0.58), ]
  up_2 = regulated_2[regulated_2$metafc > 0.58, ]
  down_2 = regulated_2[regulated_2$metafc < -0.58, ]
  names(regulated_1) <- c("GeneID", "metap", "FC1", "idx", "Entrez", "Symbol", "Description")
  names(regulated_2) <- c("GeneID", "metap", "FC2", "idx", "Entrez", "Symbol", "Description")
  combined <- merge(regulated_1, regulated_2, by = "GeneID")
  universe <- length(res_meta$GeneID) + length(res_meta2$GeneID) - length(intersect(res_meta$GeneID, res_meta2$GeneID))
  genes <- c("CD3C", "CD3D", "CD3E", "CD4", "TBX21", "IL21R", "IL1RN", "IL4I1", "IL10RA", "BMP3", "RORC")
  up1up2 = gene_overlap_test(up_1$GeneID, up_2$GeneID, universe)
  down1down2 = gene_overlap_test(down_1$GeneID, down_2$GeneID, universe)
  up1down2 = gene_overlap_test(up_1$GeneID, down_2$GeneID, universe)
  down1up2 = gene_overlap_test(down_1$GeneID, up_2$GeneID, universe)
  p2 <- ggplot(combined, aes(x=FC1, y = FC2, color = FC1, show.legend = FALSE)) +
    geom_point(size = 4, 
               shape = 16,
               show.legend = FALSE,
               alpha = 0.6) +
    bbplot::bbc_style() +
    theme(legend.position = "none") +
    scale_color_viridis_c(guide=guide_colorbar(reverse=FALSE),
                          end = 0.7,
                          begin = 0,
                          option = "inferno") + 
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    guides(
      x=guide_none(title=""),
      x.sec = guide_axis(title = paste0("Mean FC in ", disease1)),
      y = guide_none(title = ''),
      y.sec = guide_axis(title = paste0("Mean FC in ", disease2))
    ) + 
    coord_axes_inside(labels_inside = TRUE) +
    geom_label_repel(label = 
                      ifelse(combined$Symbol.x %in% genes, as.character(combined$Symbol.x), ""),
                    size=5,
                    color = "black",
                    nudge_x = 0.05,
                    box.padding = 0.5,
                    nudge_y = 0.15,
                    segment.curvature = -0.1,
                    segment.ncp = 3,
                    segment.angle = 20,
                    arrow = arrow(length = unit(0.015, "npc"))) +
    xlab(paste0("Mean FC in ", disease1)) +
    ylab(paste0("Mean FC in ", disease2)) +
    theme(legend.position = "none") +
    theme(axis.title = element_text(size = 18)) +
    annotate("text",
             x = -1.5,
             y = -1.5,
             label = sprintf("%.2e", down1down2), 
             hjust = 0, 
             vjust = 0.5, 
             colour = "black", 
             size = 7,
             fontface = "bold",
             check_overlap = TRUE) +
    annotate("text",
             x = 1.5,
             y = 1.5,
             label = sprintf("%.2e", up1up2), 
             hjust = 0, 
             vjust = 0.5, 
             colour = "black", 
             size = 7,
             fontface = "bold",
             check_overlap = TRUE) +
    annotate("text",
             x = -1.5,
             y = 1.5,
             label = sprintf("%.2e", down1up2), 
             hjust = 0, 
             vjust = 0.5, 
             colour = "black", 
             size = 7,
             fontface = "bold",
             check_overlap = TRUE) +
    annotate("text",
             x = 1.5,
             y = -1.5,
             label = sprintf("%.2e", up1down2), 
             hjust = 0, 
             vjust = 0.5, 
             colour = "black", 
             size = 7,
             fontface = "bold",
             check_overlap = TRUE)
  combined = combined[sign(combined$FC1) == sign(combined$FC2), ]
  p3 <- ggplot(enrichKEGG(combined$Entrez.x), showCategory = 10, 
               aes(-log10(p.adjust), fct_reorder(Description, -log10(p.adjust)))) + 
    geom_segment(aes(xend=0, yend = Description)) +
    geom_point(aes(color=-log10(p.adjust)), size = 10) +
    scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE),
                          end = 0.9,
                          begin = 0.4,
                          option = "inferno",
                          alpha = 0.8)  +
    bbplot::bbc_style() +
    theme(legend.position = "right", 
          legend.justification = "left",
          legend.title = element_text("Adjusted P"),
          text = element_text(size = 15)
    ) +
    xlim(c(0, 20)) +
    labs(x=bquote(~-log[10]*"(adjusted."*italic(P)*")"),
         y="",
         title="")
  title = ggdraw() + draw_label(paste0("Dataset left out: ", removed))
  intermed = cowplot::plot_grid(p2, p3)
  final_plot = cowplot::plot_grid(title, p1, intermed, nrow = 3, rel_heights = c(0.1, 1, 1))
  final_res = list(
    Meta_results_1 = res_meta,
    GSE = GSE1,
    Meta_results_2 = res_meta2,
    plot1 = p1,
    plot2 = p2,
    plot3 = p3,
    paths = shared,
    up1_up2 = up1up2,
    down1_down2 = down1down2,
    up1_down2 = up1down2,
    down1_up2 = down1up2,
    final_plot = final_plot
  )
  return(final_res)
}
