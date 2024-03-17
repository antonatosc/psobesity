library(easypackages)
libraries("tidyverse",
"ggplot2",
"clusterProfiler")
load("~/Psobesity/Psoriasis/Psoriasis.RData")
load("~/Psobesity/Obesity/Obesity.RData")
####Volcano plots
formatting <- function(result, log2_threshold, pvalue_threshold){
    result[["expression"]] <- ifelse(result[["metafc"]] > log2_threshold & -log10(result[["metap"]]) > -log10(pvalue_threshold), "Up",
        ifelse(result[["metafc"]] < -1*log2_threshold & -log10(result[["metap"]]) > -log10(pvalue_threshold) , "Down",
            ifelse(between(result[["metafc"]],-1*log2_threshold, log2_threshold) & -log10(result[["metap"]]) > -log10(pvalue_threshold), "Stable1",
                ifelse(abs(result[["metafc"]]) > log2_threshold & -log10(result[["metap"]]) < -log10(pvalue_threshold), "Stable2",
                    "Stable"))))
    return(result)
                    }
volcano_plot <- function(result, log2fc_threshold, pvalue_threshold){
    ggplot(result,
       aes(x = metafc,
           y = -log10(metap),
           color = expression)) +
  geom_point(size = 2,
             alpha = 0.7) +
  ggthemes::theme_few() +
  scale_color_manual(values=c("darkred", "black", "darkblue", "forestgreen", "darkred")) +
  geom_vline(xintercept=c(-1*log2fc_threshold, log2fc_threshold),
             col="black",
             linetype = "dashed") +
  geom_hline(yintercept=-log10(pvalue_threshold),
             col="black",
             linetype = "dashed") +
  theme(legend.position = "none") +
  xlim(c(-5, 5)) +
  labs(x=bquote("Mean"~log[2]*"(Fold Change)"),
       y=bquote("Fisher's combined -"*log[10]*"(adjusted."*italic(P)*")"),
       title="")  +
  theme(text = element_text(size = 15))
}
#####Ridgeplots
psoriasis_ridgeplot <- ridgeplot(GSE_psoriasis)
#####GSEA comparison
pathways_pso <- GSE_psoriasis@result
pathways_obs <- GSE_obesity@result
shared_pso <- pathways_pso[pathways_pso$ID %in% pathways_obs$ID, ]
shared_obs <- pathways_obs[pathways_obs$ID %in% pathways_pso$ID, ]
shared_pathways <- merge(shared_psoriasis, shared_obesity, by = "Description")
ggplot(shared_pathways, #ypo
       aes(x = NES.x,
           y = NES.y,
           color = ifelse(NES.x>0 & NES.y>0, "blue", 
                          ifelse(NES.x < 0 & NES.y < 0, "purple",
                                 ifelse(NES.x>0 & NES.y<0, "green", "red"))))) +
  geom_point(size = 4,
             alpha = 0.5) +
  stat_ellipse()
