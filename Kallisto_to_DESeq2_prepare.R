#kallisto to EdgeR via tximport

Prepare_DESeq2_12.28.17 = function(){
  
  #Load libraries, read in metadata and factor metadata ----
  library('tximport')
  library('tximportData')
  library('DESeq2')
  
  metadata = read.table("IgnaciometadataCombined.txt", sep = "\t", header = TRUE)
  metadata$path = as.character(metadata$path)
  metadata$Condition = as.factor(metadata$Condition)
  
 
  
  
  # import gene mappings (t2g) --------------------- this is a crucial section but corrently the 'mart' function fails half the time (do it by hand)
  # mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  #                          dataset = "scerevisiae_gene_ensembl",
  #                          host = 'ensembl.org')
  
  
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name", "description"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                       ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
 

 
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
  sampledata = metadata[,c("sample", "Condition","BiologicalReplicate","GrowthCondition","SequencingRun")]
  for(i in 2:length(sampledata)){
    sampledata[i] = factor(sampledata[,i])
  }
  
  # Create DESeq2 object ----
  dds <- DESeqDataSetFromTximport(txi,
                                  colData = sampledata,
                                  design = ~SequencingRun + Condition)
  # Filter based on Rows with at least 35 counts ----
  keep <- rowSums(counts(dds)) >= 35
  dds <- dds[keep,]
  
  #dds$Condition <- relevel(dds$Condition, ref = "WT Dex")
  #dds$Condition <- droplevels(dds$Condition) # only for when you are removing some of the data for the analysis
  
  # Run DESeq to create full dds object with fitted models ----
  dds = DESeq(dds)
  saveRDS(dds, "dds.rda")
  saveRDS(sampledata, "sampledata.rda")
}

