#DESeq2_Heatmaps_ect
#Will need to load in .csv files from DESeq2_Statistical and datastructure from DESeq2_prepare



DESeq2Clusterplot = function(filepath='/Users/HoltLab/Documents/R/Experiments/RNA_seq_ADH2_Nacho_12.28.17/Results/DESeq2/Sig_Dif_expressed_genes/WT_DEX_v_starve_pH_5_sig.txt', filename='WT_DEX_v_starve_pH_5_sig', numbergenes = 25){
    library("pheatmap")
    library("DESeq2")
    dds = readRDS('dds.rda') #from wherever you saved the data structure (should automatically save in prepare)----
    sampledata = readRDS('sampledata.rda')# sampledata - these both assume that these are in the current folder - can add to function w/ path if not the case
    
    # create filenames
    filename = read.table(filepath, sep = "\t", header = TRUE)
    filenameup = paste(filename,"upinstarve","_",numbergenes, sep = "") 
    filenamedown = paste(filename,"downinstarve","_",numbergenes, sep = "") 
    
    #filter SigGene lists into up and downregulated genes and cut by number desired ----
    filenamedown = filename[order(filename$log2FoldChange),][1:numbergenes,]
    
    filenameup = filename[order(filename$log2FoldChange, decreasing = TRUE),][1:numbergenes,]
   
    #create master list of all genes you selected ----
    select = c(rownames(filenameup), rownames(filenamedown))
   
    #create transformation to have a more reasonable distribution of data----
    ntd <- normTransform(dds) # can probably play with different normalizations here 
   
   
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
   
   pheatmap(cluster, cluster_rows=TRUE, show_rownames=TRUE,
            cluster_cols=TRUE, show_colnames = T, annotation_col = df)
    
}


 



