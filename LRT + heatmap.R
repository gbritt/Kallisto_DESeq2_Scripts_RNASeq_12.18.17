#DESeq2 LRT + heatmap

#when running this run 'prepare' by hand without running hte DESeq comman. Not sure if doing the LRT overwrites the Padj values from before or not. SO have a fresh start each time
dds <- DESeq(dds, test="LRT", reduced=~SequencingRun)
res <- results(dds)

res = res[order(res$padj),]
res_sig = res[res$padj < 0.05,]

res_sig_table = as.data.frame(res_sig)

ntd <- normTransform(dds)

select = rownames(res_sig_table)

cluster = assay(ntd)[select,] # making matrix with gene names and gene counts
colnames(cluster)= as.character(sampledata$sample) #giving colnames to clustered object using the 'sampledata' variable created in the prepare script

# Adding common names to ens_ids for genes ----
shared_results <- t2g[t2g[,2] %in% rownames(cluster),] # Select rownames in the 'gene mapping' matrix that are also in our clustering matrix
a = as.matrix(shared_results)
ordered_genes_mapping <- a[ order(a[,1]), ]# ordering genes mapping so that the rownames in cluster to be alphabetical
cluster <- cluster[ order(row.names(cluster)), ] # order genes in cluster to be alphabetical (will be important below)

common_and_ens_ids = matrix()  # create matrix to hold a pasted name containing common and ens_id names for a gene


for(i in 1:length(ordered_genes_mapping[,1])){
  common_and_ens_ids[i] = paste(ordered_genes_mapping[i,3], "_",ordered_genes_mapping[i,1], sep = "")
}# paste together each common name with the ens_id (will make heatmap easier to read)

rownames(cluster) = common_and_ens_ids # now our cluster rownames contain common names if available

#create annotations for top of heatmap ----
df <- as.data.frame(colData(dds)[,c("Condition","GrowthCondition")])
rownames(df) <- colnames(cluster)

pheatmap(cluster, cluster_rows=TRUE, show_rownames=F,
         cluster_cols=TRUE, show_colnames = T, annotation_col = df)