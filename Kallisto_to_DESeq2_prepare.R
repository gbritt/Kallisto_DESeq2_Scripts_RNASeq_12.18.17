#kallisto to EdgeR via tximport

Prepare_DESeq2_12.28.17 = function(){
  metadata = read.table("IgnaciometadataCombined.txt", sep = "\t", header = TRUE)
  metadata$path = as.character(metadata$path)
  metadata$Condition = as.factor(metadata$Condition)
  
  library('tximport')
  library('tximportData')
  library('DESeq2')
  
  files = metadata[,2]
  tx2gene = t2g[,1:2]
  
  
  for(i in 1:length(files)){
    files[i] = paste(files[i],sep ="","/abundance.tsv")
  } 
  
  
  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene) #, countsFromAbundance = "lengthScaledTPM" #scaled by lenght so can begin using this as your count matrix.
  
  #Apparently we really don't want to do this scaled to library size..
  
  #cts = data matrix
  # coldata = sample information
  samplesdata = metadata[,c("sample", "Condition","BiologicalReplicate","GrowthCondition","SequencingRun")]
  for(i in 2:length(samplesdata)){
    samplesdata[i] = factor(samplesdata[,i])
  }
  
  dds <- DESeqDataSetFromTximport(txi,
                                  colData = samplesdata,
                                  design = ~SequencingRun + GrowthCondition + Condition)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  #dds$Condition <- relevel(dds$Condition, ref = "WT Dex")
  #dds$Condition <- droplevels(dds$Condition) # only for when you are removing some of the data for the analysis
  
  dds = DESeq(dds)
}

