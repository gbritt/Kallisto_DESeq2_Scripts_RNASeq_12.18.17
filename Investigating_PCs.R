#PCA analysis Ingnacio RNAseq
#PC1 


x = vsd
intgroup=c("Condition")
ntop = 500


library(RColorBrewer)
library(genefilter)
library(lattice)

# pca
rv = rowVars(assay(x))
select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
pca = prcomp(t(assay(x)[select,]))
loadings <- as.data.frame(pca$rotation)

PC1_order = loadings[order(loadings[,1], decreasing = TRUE),]
PC1_list = PC1_order[,1, drop = FALSE]

PC2_order = loadings[order(loadings[,2], decreasing = TRUE),]
PC2_list = PC2_order[,2, drop = FALSE]


PC1_top = head(PC1_list,n = 20)
write.table(PC1_top,"PC1_top.txt", sep = "/t")

PC2_top = head(PC2_list,n = 20)
write.table(PC2_top,"PC2_top.txt", sep = "/t")


PC1
YLR249W   0.07878096
YNL112W   0.07705428
YGR159C   0.07696888
YDL229W   0.07387467
YAL038W   0.07060781
YDR345C   0.06984302
YDL228C   0.06501121
YOR310C   0.06476537
YLR044C   0.06364988
YAL037C-B 0.06339450

PC1_bottom = tail(PC1_list,n = 20)
write.table(PC1_bottom,"PC1_bottom.txt", sep = "/t")

PC2_bottom = tail(PC2_list,n = 20)
write.table(PC2_bottom,"PC2_bottom.txt", sep = "/t")

PC1
YGR038C-A -0.07381643
YOR348C   -0.07477779
YCR010C   -0.07504027
YLR303W   -0.07613604
YLR174W   -0.07694728
YIL057C   -0.07931134
YMR303C   -0.08077041
YJR095W   -0.08660531
YNL117W   -0.08837056
YAR035W   -0.09412765