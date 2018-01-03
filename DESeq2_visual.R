#visually exploring sample relationships, in  transform counts for computing distances or making plots
# For visual analysis we are transforming so:
#assuming you prepped the dds using Kallisto_to_DESeq2_prepare.R
vsd <- vst(dds)
#resLFC <- lfcShrink(dds, coef=2) for shrinkage in plotting heatmaps and PCA