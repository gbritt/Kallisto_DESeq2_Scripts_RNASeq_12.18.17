#kallisto to EdgeR via tximport

Prepare_DESeq2_12.28.17 = function(){
  setwd("Documents/R/Scripts/Kallisto_DESeq2_Scripts_RNASeq_12.18.17/")
  source('/Users/HoltLab/Documents/R/is.true.R')
  #Load libraries, read in metadata and factor metadata ----
  library('tximport')
  library('tximportData')
  library('DESeq2')
  
  metadata = read.table("Ignaciometadata_Individual_rep3Dex_dropped_onepH.txt", sep = "\t", header = TRUE)
  metadata$path = as.character(metadata$path)
  #metadata$Condition = as.factor(metadata$Condition)
  
 
  
  
  # import gene mappings (t2g) --------------------- this is a crucial section but corrently the 'mart' function fails half the time (do it by hand)
  
  t2g = readRDS('t2g.rda')
  #  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  #                           dataset = "scerevisiae_gene_ensembl",
  #                           host = 'ensembl.org')
  # 
  # 
  # t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
  #                                      "external_gene_name", "description"), mart = mart)
  # t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  #                      ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
 

 
  # extract file paths and append Kallisto abundance.tsv filenames -----
  files = metadata[,2]
  tx2gene = t2g[,1:2]
  
  for(i in 1:length(files)){
    files[i] = paste(files[i],sep ="","/abundance.tsv")
  } 
  
  # Use 'tximport::tximport' function to load in 'counts' data and add gene mappings to create txi counts / mappings list for use in DEseq2 ------
  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene) #, countsFromAbundance = "lengthScaledTPM" 
  #Don't use scaled by length - DESEq2 takes care of this
  
  #Create 'sampledata' variable with relevant metadata, then factor it ----
  sampledata = metadata[,c("sample_name", "GrowthCondition", "Genotype","SequencingRun")]
  for(i in 2:length(sampledata)){
    sampledata[i] = factor(sampledata[,i])
  }
  sampdata = relevel(sampledata$GrowthCondition, "Dex")
  sampdata = relevel(sampledata$Genotype, "A_WT_pH5")
  
  # Create DESeq2 object ----
  dds <- DESeqDataSetFromTximport(txi,
                                  colData = sampledata,
                                  design = ~SequencingRun + Genotype + GrowthCondition + Genotype:GrowthCondition) 
  # Filter based on Rows with at least 20 counts ----
  keep <- rowSums(counts(dds)) >= 20
  dds <- dds[keep,]
  
  #dds$Condition <- relevel(dds$Condition, ref = "WT Dex")
  #dds$Condition <- droplevels(dds$Condition) # only for when you are removing some of the data for the analysis
  
  # Run DESeq to create full dds object with fitted models ----
  dds = DESeq(dds) #it seems you can either have sequening run data or data on the WT_pH7
  saveRDS(dds, "dds_individual.rda")
  saveRDS(sampledata, "sampledata_individual.rda")
}

dds = readRDS('dds_individual.rda')
sampledata = readRDS('sampledata_individual.rda')
t2g = readRDS('t2g.rda')
#Available contrassts = 
resultsNames(dds)

# contrasts----
# the condition effect for genotype I (the main effect)
alpha = 0.05
lfcThreshold = 1

WT_stv_v_dex = results(dds, contrast=c("GrowthCondition","Stv","Dex"), alpha = alpha, lfcThreshold = lfcThreshold) #+ = up in stv vs dex for WT

    names = match(rownames(WT_stv_v_dex),t2g[,'ens_gene'])
   ordered_anno = as.data.frame(t2g[names,])

    WT_stv_v_dex$ext_gene = ordered_anno$ext_gene
    WT_stv_v_dex$description = ordered_anno$description # add gene names and descriptions
   
    WT_stv_v_dex_sig <- subset(WT_stv_v_dex, padj < alpha, log2FoldChange > lfcThreshold)
   

# the condition effect for genotype III.
# this is the main effect *plus* the interaction term
dQ_stv_v_dex = results(dds, contrast= list(c("GrowthCondition_Stv_vs_Dex","GenotypedQ_snf5.GrowthConditionStv")), alpha = alpha, lfcThreshold = lfcThreshold)

    names = match(rownames(dQ_stv_v_dex),t2g[,'ens_gene'])
    ordered_anno = as.data.frame(t2g[names,])

    dQ_stv_v_dex$ext_gene = ordered_anno$ext_gene
    dQ_stv_v_dex$description = ordered_anno$description
    
    dQ_stv_v_dex_sig <- subset(dQ_stv_v_dex, padj < alpha, log2FoldChange > lfcThreshold)

HtoA_stv_v_dex = results(dds, contrast= list(c("GrowthCondition_Stv_vs_Dex","GenotypeHtoA_snf5.GrowthConditionStv")), alpha = alpha, lfcThreshold = lfcThreshold)

    names = match(rownames(HtoA_stv_v_dex),t2g[,'ens_gene'])
    ordered_anno = as.data.frame(t2g[names,])

    HtoA_stv_v_dex$ext_gene = ordered_anno$ext_gene
    HtoA_stv_v_dex$description = ordered_anno$description
    
    HtoA_stv_v_dex_sig <- subset(HtoA_stv_v_dex, padj < alpha, log2FoldChange > lfcThreshold)

Null_stv_v_dex = results(dds, contrast= list(c("GrowthCondition_Stv_vs_Dex","GenotypeNull.GrowthConditionStv")), alpha = alpha, lfcThreshold = lfcThreshold)

    names = match(rownames(Null_stv_v_dex),t2g[,'ens_gene'])
    ordered_anno = as.data.frame(t2g[names,])

    Null_stv_v_dex$ext_gene = ordered_anno$ext_gene
    Null_stv_v_dex$description = ordered_anno$description
    
    Null_stv_v_dex_sig <- subset(Null_stv_v_dex, padj < alpha, log2FoldChange > lfcThreshold)

#---- turning DE genes into Log2FC values for plotting - use with other scripts---- 
WT_stv_v_dex_df = as.data.frame(WT_stv_v_dex)
dQ_stv_v_dex_df = as.data.frame(dQ_stv_v_dex)
HtoA_stv_v_dex_df = as.data.frame(HtoA_stv_v_dex)
Null_stv_v_dex_df = as.data.frame(Null_stv_v_dex)

a = as.data.frame(WT_stv_v_dex_df$log2FoldChange)
rownames(a) = paste(rownames(WT_stv_v_dex_df), '_', WT_stv_v_dex_df$ext_gene, sep = "")
colnames(a) = 'WT'

b =as.data.frame(dQ_stv_v_dex_df$log2FoldChange)
rownames(b) = paste(rownames(dQ_stv_v_dex_df), '_', WT_stv_v_dex_df$ext_gene, sep = "")
colnames(b) = 'dQ'

c =as.data.frame(HtoA_stv_v_dex_df$log2FoldChange)
rownames(c) = paste(rownames(HtoA_stv_v_dex_df), '_', WT_stv_v_dex_df$ext_gene, sep = "")
colnames(c) = 'HtoA'

d =as.data.frame(Null_stv_v_dex_df$log2FoldChange)
rownames(d) = paste(rownames(Null_stv_v_dex_df), '_', WT_stv_v_dex_df$ext_gene, sep = "")
colnames(d) = 'Null'

DE_Genes_in_Starvation = cbind(a,b,c,d)
saveRDS(DE_Genes_in_Starvation, "DE_Genes_in_Starvation_master_uncut.rda")

#end----

## the interaction term for condition effect in genotype III vs genotype I.
# this tests if the condition effect is different in III compared to I
dQ_WT_interaction_stv_v_dex = results(dds, name="GenotypedQ_snf5.GrowthConditionStv", alpha = alpha, lfcThreshold = lfcThreshold)

    names = match(rownames(dQ_WT_interaction_stv_v_dex),t2g[,'ens_gene'])
    ordered_anno = as.data.frame(t2g[names,])
    
    dQ_WT_interaction_stv_v_dex$ext_gene = ordered_anno$ext_gene
    dQ_WT_interaction_stv_v_dex$description = ordered_anno$description
    
    dQ_WT_interaction_stv_v_dex_sig <- subset(dQ_WT_interaction_stv_v_dex, padj < alpha, log2FoldChange > lfcThreshold)

HtoA_WT_interaction_stv_v_dex = results(dds, name="GenotypeHtoA_snf5.GrowthConditionStv", alpha = alpha, lfcThreshold = lfcThreshold)

    names = match(rownames(HtoA_WT_interaction_stv_v_dex),t2g[,'ens_gene'])
    ordered_anno = as.data.frame(t2g[names,])
    
    HtoA_WT_interaction_stv_v_dex$ext_gene = ordered_anno$ext_gene
    HtoA_WT_interaction_stv_v_dex$description = ordered_anno$description
    
    HtoA_WT_interaction_stv_v_dex_sig <- subset(HtoA_WT_interaction_stv_v_dex, padj < alpha, log2FoldChange > lfcThreshold)

Null_WT_interaction_stv_v_dex = results(dds, name="GenotypeNull.GrowthConditionStv", alpha = alpha, lfcThreshold = lfcThreshold)

    names = match(rownames(Null_WT_interaction_stv_v_dex),t2g[,'ens_gene'])
    ordered_anno = as.data.frame(t2g[names,])
    
    Null_WT_interaction_stv_v_dex$ext_gene = ordered_anno$ext_gene
    Null_WT_interaction_stv_v_dex$description = ordered_anno$description
    
    Null_WT_interaction_stv_v_dex_sig <- subset(Null_WT_interaction_stv_v_dex, padj < alpha, log2FoldChange > lfcThreshold)


## the interaction term for condition effect in genotype III vs genotype II.
# this tests if the condition effect is different in III compared to II
dQ_HtoA_interaction_stv_v_dex = results(dds, contrast= list(c("GenotypeHtoA_snf5.GrowthConditionStv","Genotype_dQ_snf5_vs_A_WT_pH5")), alpha = alpha, lfcThreshold = lfcThreshold)

    names = match(rownames(dQ_HtoA_interaction_stv_v_dex),t2g[,'ens_gene'])
    ordered_anno = as.data.frame(t2g[names,])
    
    dQ_HtoA_interaction_stv_v_dex$ext_gene = ordered_anno$ext_gene
    dQ_HtoA_interaction_stv_v_dex$description = ordered_anno$description
    
    dQ_HtoA_interaction_stv_v_dex_sig <- subset(dQ_HtoA_interaction_stv_v_dex, padj < alpha, log2FoldChange > lfcThreshold)



#filtering and saving values ----
results_path = '/Users/HoltLab/Documents/R/Experiments/RNA_seq_ADH2_Nacho_12.28.17/Results/DESeq2/DE_genes_individual_comparisons/'
# Writing the list to CSV ----

write.table(as.data.frame(WT_stv_v_dex_sig), file = paste(results_path,"WT_stv_v_dex.txt", sep = ""), sep ="\t")
write.table(as.data.frame(dQ_stv_v_dex_sig), file = paste(results_path,"dQ_stv_v_dex.txt", sep = ""), sep ="\t")
write.table(as.data.frame(HtoA_stv_v_dex_sig), file = paste(results_path,"HtoA_stv_v_dex.txt", sep = ""), sep ="\t")
write.table(as.data.frame(Null_stv_v_dex_sig), file = paste(results_path,"Null_stv_v_dex.txt", sep = ""), sep ="\t")
write.table(as.data.frame(dQ_WT_interaction_stv_v_dex_sig), file = paste(results_path,"dQ_WT_interaction_stv_v_dex.txt", sep = ""), sep ="\t")
write.table(as.data.frame(HtoA_WT_interaction_stv_v_dex_sig), file = paste(results_path,"HtoA_WT_interaction_stv_v_dex.txt", sep = ""), sep ="\t")
write.table(as.data.frame(Null_WT_interaction_stv_v_dex_sig), file = paste(results_path,"Null_WT_interaction_stv_v_dex.txt", sep = ""), sep ="\t")
write.table(as.data.frame(dQ_HtoA_interaction_stv_v_dex_sig), file = paste(results_path,"dQ_HtoA_interaction_stv_v_dex.txt", sep = ""), sep ="\t")
# want to read it back in into matrices? Need to decide what to do next. Use genes for heatmap? Venn Diagram? GO enrichment? Just list them?
