#DESeq2_Heatmaps_ect
#Will need to load in .csv files from DESeq2_Statistical 

library("pheatmap")
DESeq2Clusterplot = function(filepath='/Users/gbrittingham/Documents/R/R_from macpro/Kallisto_DESeq2_Scripts_RNASeq_12.18.17/Results/DESeq2/Sig_Dif_expressed_genes/WT_DEX_v_starve_pH_5_sig.txt', filename='WT_DEX_v_starve_pH_5_sig', numbergenes = 25){
    load('dds.rda') #from wherever you saved the data structure
    filename = read.table(filepath, sep = "\t", header = TRUE)
    filenameup = paste(filename,"upinstarve","_",numbergenes, sep = "") 
    filenamedown = paste(filename,"downinstarve","_",numbergenes, sep = "") 
    
    filenamedown = filename[order(filename$log2FoldChange),][1:numbergenes,]
    
   filenameup = filename[order(filename$log2FoldChange, decreasing = TRUE),][1:numbergenes,]
   
   select = c(rownames(filenameup), rownames(filenamedown))
   
   #create transformation
   ntd <- normTransform(dds)
   
   
   cluster = assay(ntd)[select,] # making an object out of genes to be clusterd
   colnames(cluster)= as.character(samples$sample) #giving colnames to clustered object
   
   
   shared_results <- t2g[t2g[,2] %in% rownames(cluster),] # subsetting genes to transcripts mapping to add gene names to heatmap
   a = as.matrix(shared_results)
   ordered_genes_mapping <- a[ order(a[,1]), ]# ordering genes mapping
   
   cluster <- cluster[ order(row.names(cluster)), ] # orderign cluster rownames so that they are in     the same order as the gene names we have selected
   common_and_ens_ids = matrix()
   for(i in 1:length(ordered_genes_mapping[,1])){
       common_and_ens_ids[i] = paste(ordered_genes_mapping[i,3], "_",ordered_genes_mapping[i,1], sep = "")
   }
   #pasting the common gene name to the ens_id of each gene so we can have both for the heatmap
   
   rownames(cluster) = common_and_ens_ids
   
   df <- as.data.frame(colData(dds)[,c("Condition","GrowthCondition")])
   rownames(df) <- colnames(cluster)
   
    
}
WT_DEX_v_starve_pH_5_sig = read.table("WT_DEX_v_starve_pH_5_sig.txt", sep = "\t", header = TRUE)
WT_Dex_v_Starve_pH_7_sig = read.table("WT_Dex_v_Starve_pH_7_sig.txt", sep = "\t", header = TRUE)

WT_DEX_v_starve_pH_5_sig_upinstarve_25 = WT_DEX_v_starve_pH_5_sig[order(WT_DEX_v_starve_pH_5_sig$log2FoldChange),][1:25,]
WT_DEX_v_starve_pH_5_sig_downinstarve_25 = WT_DEX_v_starve_pH_5_sig[order(WT_DEX_v_starve_pH_5_sig$log2FoldChange, decreasing = TRUE),][1:25,]

select = c(rownames(WT_DEX_v_starve_pH_5_sig_upinstarve_25), rownames(WT_DEX_v_starve_pH_5_sig_downinstarve_25))

#create transformation
ntd <- normTransform(dds)
df <- as.data.frame(colData(dds)[,"Condition"])
#df <- as.data.frame(colData(dds)[,c("Condition","GrowthCondition")])

cluster = assay(ntd)[select,] # making an object out of genes to be clusterd
colnames(cluster)= as.character(samples$sample) #giving colnames to clustered object


shared_results <- t2g[t2g[,2] %in% rownames(cluster),] # subsetting genes to transcripts mapping to add gene names to heatmap
a = as.matrix(shared_results)
ordered_genes_mapping <- a[ order(a[,1]), ]# ordering genes mapping

cluster <- cluster[ order(row.names(cluster)), ] # orderign cluster rownames so that they are in the same order as the gene names we have selected
common_and_ens_ids = matrix()
for(i in 1:length(ordered_genes_mapping[,1])){
  common_and_ens_ids[i] = paste(ordered_genes_mapping[i,3], "_",ordered_genes_mapping[i,1], sep = "")
}
#pasting the common gene name to the ens_id of each gene so we can have both for the heatmap

rownames(cluster) = common_and_ens_ids

df <- as.data.frame(colData(dds)[,c("Condition","GrowthCondition")])
rownames(df) <- colnames(cluster)

pheatmap(cluster, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, show_colnames = T, annotation_col = df) #annotation_col=df,



