library(easypackages)
libraries("tidyverse",
"ggplot2",
"clusterProfiler",
'data.table',
'ggrepel',
'ggh4x',
'clusterProfiler')
load("~/Psobesity/Psoriasis/Psoriasis.RData")
load("~/Psobesity/Obesity/Obesity.RData")
regulated_psoriasis <- res_psoriasis[(res_psoriasis$metap < 0.01) & !between(res_psoriasis$metafc, -0.58, 0.58), ]
regulated_Obesity <- res_Obesity[(res_Obesity$metap < 0.01) & !between(res_Obesity$metafc, -0.58, 0.58), ]
names(regulated_Obesity) <- c("GeneID", "metap", "FC in SAT", "idx", "Entrez", "Symbol", "Description")
names(regulated_psoriasis) <- c("GeneID", "metap", "FC in Skin", "idx", "Entrez", "Symbol", "Description")
combined <- merge(regulated_psoriasis, regulated_Obesity, by = "GeneID")
combined <- combined[combined$`FC in SAT`*combined$`FC in Skin` > 0, ]
genes <- c("CD3C", "CD3D", "CD3E", "CD4", "TBX21", "IL21R", "IL1RN", "IL4I1", "IL10RA", "BMP3", "RORC")
p1 <- ggplot(combined, aes(x=`FC in SAT`, y = `FC in Skin`, color = `FC in Skin`, , show.legend = FALSE)) +
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
    x.sec = guide_axis(title = 'Mean FC in SAT'),
    y = guide_none(title = ''),
    y.sec = guide_axis(title = "Mean FC in Skin")
  ) + 
  coord_axes_inside(labels_inside = TRUE) +
  geom_text_repel(label = 
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
  xlab("Mean FC in SAT") +
  ylab("Mean FC in Skin") +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 18)) + 
  annotate("text",
           x = -0.5,
           y = -0.25,
           label = bquote(~n[genes]*"=49, P-value=7.1*"*10^-7), 
           hjust = 0, 
           vjust = 0.5, 
           colour = "black", 
           size = 5,
           check_overlap = TRUE) +
  annotate("text",
           x = 1.25,
           y = 2.5,
           label = bquote(~n[genes]*"=170, P=6.07*"*10^-65), 
           hjust = 0, 
           vjust = 0.5, 
           colour = "black", 
           size = 5,
           check_overlap = TRUE)
p2 <- ggplot(enrichKEGG(combined$Entrez.x), showCategory = 10, 
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