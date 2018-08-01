#DESeq2 differential gene expression analysis
#Assuming you ran Kallisto_to_DESeq2_prepare.R
# Currently not a function
library("IHW")
library("DESeq2")


dds = readRDS('dds.rda')
sampledata = readRDS('sampledata.rda')
t2g = readRDS('t2g.rda')

alpha = 0.05
lfcThreshold = 1



# Create Wald Test contrasts ----


WT_Starve_v_Dex_pH_7 = results(dds, contrast=c("Condition", "WT_Stv_pH7","WT_Dex"), alpha = alpha, lfcThreshold = lfcThreshold) #+ = up in stv vs dex for WT

    names = match(rownames(WT_Starve_v_Dex_pH_7),t2g[,'ens_gene'])
    ordered_anno = as.data.frame(t2g[names,])
    
    WT_Starve_v_Dex_pH_7$ext_gene = ordered_anno$ext_gene
    WT_Starve_v_Dex_pH_7$description = ordered_anno$description # add gene names and descriptions
    
    WT_Starve_v_Dex_pH_7_sig <- subset(WT_Starve_v_Dex_pH_7, padj < alpha, log2FoldChange > lfcThreshold)
    
#
WT_Starve_v_Dex = results(dds, contrast=c("Condition","WT_Stv_pH5","WT_Dex") , alpha = alpha, lfcThreshold = lfcThreshold)

    names = match(rownames(WT_Starve_v_Dex),t2g[,'ens_gene'])
    ordered_anno = as.data.frame(t2g[names,])
    
    WT_Starve_v_Dex$ext_gene = ordered_anno$ext_gene
    WT_Starve_v_Dex$description = ordered_anno$description # add gene names and descriptions
    
    WT_Starve_v_Dex_sig <- subset(WT_Starve_v_Dex, padj < alpha, log2FoldChange > lfcThreshold)

#
dQ_Starve_v_Dex = results(dds, contrast=c("Condition","dQ_snf5_Stv_pH5", "dQ_snf5_Dex"), alpha = alpha, lfcThreshold = lfcThreshold)
    names = match(rownames(dQ_Starve_v_Dex),t2g[,'ens_gene'])
    ordered_anno = as.data.frame(t2g[names,])
    
    dQ_Starve_v_Dex$ext_gene = ordered_anno$ext_gene
    dQ_Starve_v_Dex$description = ordered_anno$description # add gene names and descriptions
    
    dQ_Starve_v_Dex_sig <- subset(dQ_Starve_v_Dex, padj < alpha, log2FoldChange > lfcThreshold)

Null_Starve_v_Dex = results(dds, contrast=c("Condition","Null_Stv_pH5","Null_Dex"), alpha = alpha, lfcThreshold = lfcThreshold)

  names = match(rownames(Null_Starve_v_Dex),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  Null_Starve_v_Dex$ext_gene = ordered_anno$ext_gene
  Null_Starve_v_Dex$description = ordered_anno$description # add gene names and descriptions
  
  Null_Starve_v_Dex_sig <- subset(Null_Starve_v_Dex, padj < alpha, log2FoldChange > lfcThreshold)


HtoA_Starve_v_Dex = results(dds, contrast=c("Condition","HtoA_snf5_Stv_pH5", "HtoA_snf5_Dex"), alpha = alpha, lfcThreshold = lfcThreshold)

  names = match(rownames(HtoA_Starve_v_Dex),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  HtoA_Starve_v_Dex$ext_gene = ordered_anno$ext_gene
  HtoA_Starve_v_Dex$description = ordered_anno$description # add gene names and descriptions
  
  HtoA_Starve_v_Dex_sig <- subset(HtoA_Starve_v_Dex, padj < alpha, log2FoldChange > lfcThreshold)
  


#Results path
results_path = '/Users/HoltLab/Documents/R/Experiments/RNA_seq_ADH2_Nacho_12.28.17/Results/DESeq2/DE_genes_combined/'
# Writing the list to CSV ----
write.table(as.data.frame(WT_Starve_v_Dex_pH_7_sig), file = paste(results_path,"WT_Starve_v_Dex_pH_7.txt", sep = ""), sep ="\t")
write.table(as.data.frame(WT_Starve_v_Dex_sig), file= paste(results_path, "WT_Starve_v_Dex.txt", sep = ""), sep ="\t")
write.table(as.data.frame(dQ_Starve_v_Dex_sig), file= paste(results_path, "dQ_Starve_v_Dex.txt", sep = ""), sep ="\t")
write.table(as.data.frame(Null_Starve_v_Dex_sig), file= paste(results_path, "Null_Starve_v_Dex.txt", sep = ""), sep ="\t")
write.table(as.data.frame(HtoA_Starve_v_Dex_sig), file= paste(results_path, "HtoA_Starve_v_Dex.txt", sep = ""), sep ="\t")

# Without filtering and with lfcShrinking (needed for volcano plot) ----

WT_Starve_v_Dex_pH_7 = lfcShrink(dds, contrast=c("Condition", "WT_Stv_pH7","WT_Dex"), cooksCutoff=FALSE, independentFiltering=FALSE) #+ = up in stv vs dex for WT

names = match(rownames(WT_Starve_v_Dex_pH_7),t2g[,'ens_gene'])
ordered_anno = as.data.frame(t2g[names,])
  
  WT_Starve_v_Dex_pH_7$ext_gene = ordered_anno$ext_gene
  WT_Starve_v_Dex_pH_7$description = ordered_anno$description # add gene names and descriptions
  
  WT_Starve_v_Dex_pH_7_Volcano <- subset(WT_Starve_v_Dex_pH_7, padj < alpha, log2FoldChange > lfcThreshold)

#
WT_Starve_v_Dex = lfcShrink(dds, contrast=c("Condition","WT_Stv_pH5","WT_Dex") , cooksCutoff=FALSE, independentFiltering=FALSE)

  names = match(rownames(WT_Starve_v_Dex),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  WT_Starve_v_Dex$ext_gene = ordered_anno$ext_gene
  WT_Starve_v_Dex$description = ordered_anno$description # add gene names and descriptions
  
  WT_Starve_v_Dex_Volcano <- subset(WT_Starve_v_Dex, padj < alpha, log2FoldChange > lfcThreshold)

#
dQ_Starve_v_Dex = lfcShrink(dds, contrast=c("Condition","dQ_snf5_Stv_pH5", "dQ_snf5_Dex"), cooksCutoff=FALSE, independentFiltering=FALSE)
  names = match(rownames(dQ_Starve_v_Dex),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  dQ_Starve_v_Dex$ext_gene = ordered_anno$ext_gene
  dQ_Starve_v_Dex$description = ordered_anno$description # add gene names and descriptions
  
  dQ_Starve_v_Dex_Volcano <- subset(dQ_Starve_v_Dex, padj < alpha, log2FoldChange > lfcThreshold)

Null_Starve_v_Dex = lfcShrink(dds, contrast=c("Condition","Null_Stv_pH5","Null_Dex"), cooksCutoff=FALSE, independentFiltering=FALSE)

  names = match(rownames(Null_Starve_v_Dex),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  Null_Starve_v_Dex$ext_gene = ordered_anno$ext_gene
  Null_Starve_v_Dex$description = ordered_anno$description # add gene names and descriptions
  
  Null_Starve_v_Dex_Volcano <- subset(Null_Starve_v_Dex, padj < alpha, log2FoldChange > lfcThreshold)


HtoA_Starve_v_Dex = lfcShrink(dds, contrast=c("Condition","HtoA_snf5_Stv_pH5", "HtoA_snf5_Dex"), cooksCutoff=FALSE, independentFiltering=FALSE)

  names = match(rownames(HtoA_Starve_v_Dex),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  HtoA_Starve_v_Dex$ext_gene = ordered_anno$ext_gene
  HtoA_Starve_v_Dex$description = ordered_anno$description # add gene names and descriptions
  
  HtoA_Starve_v_Dex_Volcano <- subset(HtoA_Starve_v_Dex, padj < alpha, log2FoldChange > lfcThreshold)
  
 #comparing within starvation ####----- 
dQ_v_WT_Starve = lfcShrink(dds, contrast=c("Condition","WT_Stv_pH5", "dQ_snf5_Stv_pH5"), cooksCutoff=FALSE, independentFiltering=FALSE)
  
  names = match(rownames(dQ_v_WT_Starve),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  dQ_v_WT_Starve$ext_gene = ordered_anno$ext_gene
  dQ_v_WT_Starve$description = ordered_anno$description # add gene names and descriptions
  
  dQ_v_WT_Starve_sig <- subset(dQ_v_WT_Starve, padj < alpha, log2FoldChange > lfcThreshold)
  
dQ_v_WT_Starve_pH7 = lfcShrink(dds, contrast=c("Condition","WT_Stv_pH7", "dQ_snf5_Stv_pH5"), cooksCutoff=FALSE, independentFiltering=FALSE)
  
  names = match(rownames(dQ_v_WT_Starve_pH7),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  dQ_v_WT_Starve_pH7$ext_gene = ordered_anno$ext_gene
  dQ_v_WT_Starve_pH7$description = ordered_anno$description # add gene names and descriptions
  
  dQ_v_WT_Starve_pH7_sig <- subset(dQ_v_WT_Starve_pH7, padj < alpha, log2FoldChange > lfcThreshold)
  
dQ_v_HtoA_Starve = lfcShrink(dds, contrast=c("Condition","HtoA_snf5_Stv_pH5", "dQ_snf5_Stv_pH5"), cooksCutoff=FALSE, independentFiltering=FALSE)
  
  names = match(rownames(dQ_v_HtoA_Starve),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  dQ_v_HtoA_Starve$ext_gene = ordered_anno$ext_gene
  dQ_v_HtoA_Starve$description = ordered_anno$description # add gene names and descriptions
  
  dQ_v_WT_Starve_pH7_sig <- subset(dQ_v_HtoA_Starve, padj < alpha, log2FoldChange > lfcThreshold)
  
dQ_v_Null_Starve = lfcShrink(dds, contrast=c("Condition","Null_Stv_pH5", "dQ_snf5_Stv_pH5"), cooksCutoff=FALSE, independentFiltering=FALSE)
  
  names = match(rownames(dQ_v_Null_Starve),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  dQ_v_Null_Starve$ext_gene = ordered_anno$ext_gene
  dQ_v_Null_Starve$description = ordered_anno$description # add gene names and descriptions
  
  dQ_v_WT_Starve_pH7_sig <- subset(dQ_v_Null_Starve, padj < alpha, log2FoldChange > lfcThreshold)
  
WT_v_Null_Starve = lfcShrink(dds, contrast=c("Condition","WT_Stv_pH5", "Null_Stv_pH5"), cooksCutoff=FALSE, independentFiltering=FALSE)
  
  names = match(rownames(WT_v_Null_Starve),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  WT_v_Null_Starve$ext_gene = ordered_anno$ext_gene
  WT_v_Null_Starve$description = ordered_anno$description # add gene names and descriptions
  
  dQ_v_WT_Starve_pH7_sig <- subset(WT_v_Null_Starve, padj < alpha, log2FoldChange > lfcThreshold)
  
WT_pH4_v_7_Starve = lfcShrink(dds, contrast=c("Condition","WT_Stv_pH5","WT_Stv_pH7") , cooksCutoff=FALSE, independentFiltering=FALSE)
  
  names = match(rownames(WT_pH4_v_7_Starve),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  WT_pH4_v_7_Starve$ext_gene = ordered_anno$ext_gene
  WT_pH4_v_7_Starve$description = ordered_anno$description # add gene names and descriptions
  
  WT_pH4_v_7_Starve_Volcano <- subset(WT_pH4_v_7_Starve, padj < alpha, log2FoldChange > lfcThreshold) 

# ---- Create DE genes log2FC list for use with Heatmap ----
  
  
  #need to add the following from individual!!
  #---- turning DE genes into Log2FC values for plotting - use with other scripts---- 
  WT_Starve_v_Dex_df = as.data.frame(WT_Starve_v_Dex)
  dQ_Starve_v_Dex_df = as.data.frame(dQ_Starve_v_Dex)
  HtoA_Starve_v_Dex_df = as.data.frame(HtoA_Starve_v_Dex)
  Null_Starve_v_Dex_df = as.data.frame(Null_Starve_v_Dex)
  WT_Starve_v_Dex_pH_7_df = as.data.frame(WT_Starve_v_Dex_pH_7)
  
  a = as.data.frame(WT_Starve_v_Dex_df$log2FoldChange)
  rownames(a) = paste(rownames(WT_Starve_v_Dex_df), '_', WT_Starve_v_Dex_df$ext_gene, sep = "")
  colnames(a) = 'WT'
  
  b =as.data.frame(dQ_Starve_v_Dex_df$log2FoldChange)
  rownames(b) = paste(rownames(dQ_Starve_v_Dex_df), '_', dQ_Starve_v_Dex_df$ext_gene, sep = "")
  colnames(b) = 'dQ'
  
  c =as.data.frame(HtoA_Starve_v_Dex_df$log2FoldChange)
  rownames(c) = paste(rownames(HtoA_Starve_v_Dex_df), '_', HtoA_Starve_v_Dex_df$ext_gene, sep = "")
  colnames(c) = 'HtoA'
  
  d =as.data.frame( Null_Starve_v_Dex_df$log2FoldChange)
  rownames(d) = paste(rownames( Null_Starve_v_Dex_df), '_',  Null_Starve_v_Dex_df$ext_gene, sep = "")
  colnames(d) = 'Null'
  
  e =as.data.frame(  WT_Starve_v_Dex_pH_7_df$log2FoldChange)
  rownames(e) = paste(rownames(  WT_Starve_v_Dex_pH_7_df), '_',   WT_Starve_v_Dex_pH_7_df$ext_gene, sep = "")
  colnames(e) = 'WT_pH7'
  
  DE_Genes_in_Starvation = cbind(a,b,c,d,e)
  saveRDS(DE_Genes_in_Starvation, "DE_Genes_in_Starvation_combined_master.rda")
  
  
  # Without the '_' and common gene name ----
  
  WT_Starve_v_Dex_df = as.data.frame(WT_Starve_v_Dex)
  dQ_Starve_v_Dex_df = as.data.frame(dQ_Starve_v_Dex)
  HtoA_Starve_v_Dex_df = as.data.frame(HtoA_Starve_v_Dex)
  Null_Starve_v_Dex_df = as.data.frame(Null_Starve_v_Dex)
  WT_Starve_v_Dex_pH_7_df = as.data.frame(WT_Starve_v_Dex_pH_7)
  
  a = as.data.frame(WT_Starve_v_Dex_df$log2FoldChange)
  rownames(a) = rownames(WT_Starve_v_Dex_df)
  colnames(a) = 'WT'
  
  b =as.data.frame(dQ_Starve_v_Dex_df$log2FoldChange)
  rownames(b) = rownames(dQ_Starve_v_Dex_df)
  colnames(b) = 'dQ'
  
  c =as.data.frame(HtoA_Starve_v_Dex_df$log2FoldChange)
  rownames(c) = rownames(HtoA_Starve_v_Dex_df)
  colnames(c) = 'HtoA'
  
  d =as.data.frame( Null_Starve_v_Dex_df$log2FoldChange)
  rownames(d) = rownames( Null_Starve_v_Dex_df)
  colnames(d) = 'Null'
  
  e =as.data.frame(  WT_Starve_v_Dex_pH_7_df$log2FoldChange)
  rownames(e) = rownames(  WT_Starve_v_Dex_pH_7_df)
  colnames(e) = 'WT_pH7'
  
  DE_Genes_in_Starvation = cbind(a,b,c,d,e)
  saveRDS(DE_Genes_in_Starvation, "DE_Genes_in_Starvation_combined_master.rda")
  
  gn.most.sign <- rownames(res.DESeq2)[1]
  gn.most.diff.val <- counts(dds.norm, normalized=T)[gn.most.sign,]
  barplot(gn.most.diff.val, col=col.pheno.selected, main=gn.most.sign, las=2, cex.names=0.5)
  
  res_table = counts(dds, normalized=T) #gene counts results table
  colnames(res_table) =  dds$sample
  res_ordered = res_table[,c(1,10,23,2,11,24,3,12,25,4,13,26,5,14,19,27,6,15,20,28,7,16,21,29,8,17,22,30,9,18,31)]
  #for when I dropped a bunch of samples
  
  res_ordered = res_table[,c(1,7,19,2,8,20,3,9,21,4,10,22,5,11,23,6,15,27,12,16,24,13,17,25,14,18,26)] #for when I dropped a bunch of samples

  
  DEstarveGenes <- WT_Starve_v_Dex_sig[order(WT_Starve_v_Dex_sig$pvalue),] # can use these with heatmap

  dQ_genes_expr = dQ_Starve_v_Dex_sig[order(dQ_Starve_v_Dex_sig$pvalue),]
  
  
  
  sum(a[,1])
  sum(a[,2])
  sum(a[,3])
  sum(a[,4])
  sum(a[,5])
  
  
  fold_induction = read.table("Fold_induction_starvation_response_genes.txt")
  
  a = (fold_induction > 5)
  WT_5_names = rownames(fold_induction[a[,1],])
  WT_7_names = rownames(fold_induction[a[,2],])
  dQ_names = rownames(fold_induction[a[,3],])
  HtoA_names = rownames(fold_induction[a[,4],])
  Null_names = rownames(fold_induction[a[,5],])

write.table(WT_5_names, "WT_5_names.txt", sep = "\t")
write.table(WT_7_names, "WT_7_names.txt", sep = "\t")
write.table(dQ_names, "dQ_Names.txt", sep = "\t")
write.table(HtoA_names,"HtoA_names.txt", sep = "\t")
write.table(Null_names, "Null_names.txt", sep = "\t")


fold_induction_ratio = cbind(fold_induction[,2]/fold_induction[,1], fold_induction[,3]/fold_induction[,1], fold_induction[,4]/fold_induction[,1], fold_induction[,5]/fold_induction[,1])

fold_induction_ratio[is.na(fold_induction_ratio)] = 0

fold_induction_ratio[(fold_induction_ratio) == Inf] = 0
colnames(fold_induction_ratio) =  colnames(fold_induction[,2:5])
rownames(fold_induction_ratio) = rownames(fold_induction)

subset_induction = fold_induction_ratio > 0.85

WT_pH5_cutoff = rownames(subset_induction)
pH7_cutoff = rownames(subset_induction)[subset_induction[,1]]
dQ_cutoff = rownames(subset_induction)[subset_induction[,2]]
HtoA_cutoff = rownames(subset_induction)[subset_induction[,3]]
Null_cutoff = rownames(subset_induction)[subset_induction[,4]]

89
19
5
11
23

v_left = read.table('/Users/HoltLab/Dropbox/ADH2_Switch_snf_RNAseq_12.28.17/volcano_leftside.txt', sep = "\t")
v_right = read.table('/Users/HoltLab/Dropbox/ADH2_Switch_snf_RNAseq_12.28.17/volcano_rightside.txt', sep = "\t")

v_left = as.character(v_left[,1])
v_right = as.character(v_right[,1])

t2g_sub_left = t2g$ext_gene %in% v_left
t2g_sub_right = t2g$ext_gene %in% v_right

left_genes = t2g[t2g_sub_left,]
right_genes = t2g[t2g_sub_right,]

write.table(left_genes, '/Users/HoltLab/Dropbox/ADH2_Switch_snf_RNAseq_12.28.17/volcano_annot_left.txt', sep = "\t")
write.table(right_genes, '/Users/HoltLab/Dropbox/ADH2_Switch_snf_RNAseq_12.28.17/volcano_annot_right.txt', sep = "\t")