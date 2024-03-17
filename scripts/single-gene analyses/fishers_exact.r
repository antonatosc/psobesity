gene_overlap_test = function(gene_list1, gene_list2, background){
  A = length(gene_list1)
  B = length(gene_list2)
  C = background
  AB = length(intersect(gene_list1, gene_list2))
  pval = phyper(AB-1, B, C-B, A, lower.tail = FALSE)
  writeLines(paste0("Length of list A: ", A, ".\nLength of list B: ", B, ".\n", AB, " matching from lists A and B.\nP-value:",pval)) #typo
}
