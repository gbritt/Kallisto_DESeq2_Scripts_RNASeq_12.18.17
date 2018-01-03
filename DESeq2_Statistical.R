#DESeq2 differential gene expression analysis

#assuming you ran Kallisto_to_DESeq2_prepare.R
WT_Dex_v_Starve_pH_7 = results(dds, contrast=c("Condition","WT_Dex", "WT_Stv_pH7"))
WT_Dex_v_Starve_pH_5 = results(dds, contrast=c("Condition","WT_Dex", "WT_Stv_pH5"))
dQ_snf5_Dex_v_Starve_pH_5 = results(dds, contrast=c("Condition","dQ_snf5_Dex", "dQ_snf5_Stv_pH5" ))
Null_Dex_v_Starve = results(dds, contrast=c("Condition","Null_Dex", "Null_Stv_pH5"))
HtoA_snf5_Dex_v_Starve_pH_5 = results(dds, contrast=c("Condition","HtoA_snf5_Dex", "HtoA_snf5_Stv_pH5"))


#now ordering them from greatest to least FC difference
WT_Dex_v_Starve_pH_7 = WT_Dex_v_Starve_pH_7[order(WT_Dex_v_Starve_pH_7$padj),]
WT_Dex_v_Starve_pH_5 = WT_Dex_v_Starve_pH_5[order(WT_Dex_v_Starve_pH_5$padj),]
dQ_snf5_Dex_v_Starve_pH_5 = dQ_snf5_Dex_v_Starve_pH_5[order(dQ_snf5_Dex_v_Starve_pH_5$padj),]
Null_Dex_v_Starve = Null_Dex_v_Starve[order(Null_Dex_v_Starve$padj),]
HtoA_snf5_Dex_v_Starve_pH_5 = HtoA_snf5_Dex_v_Starve_pH_5[order(HtoA_snf5_Dex_v_Starve_pH_5$padj),]
#Can also use interactions: Interaction terms can be added to the design formula, in order to test, for example, 
#if the log2 fold change attributable to a given condition is different based on another 
#factor, for example if the condition effect differs across genotype.

# Now filtering by padj value
  WT_Dex_v_Starve_pH_7_sig = WT_Dex_v_Starve_pH_7[WT_Dex_v_Starve_pH_7$padj < 0.05,]
  WT_DEX_v_starve_pH_5_sig = WT_Dex_v_Starve_pH_5[WT_Dex_v_Starve_pH_5$padj < 0.05,]
  dQ_snf5_Dex_v_Starve_pH_5_sig = dQ_snf5_Dex_v_Starve_pH_5[dQ_snf5_Dex_v_Starve_pH_5$padj < 0.05,]
  Null_Dex_v_Starve_sig = Null_Dex_v_Starve[Null_Dex_v_Starve$padj < 0.05,]
  HtoA_snf5_Dex_v_Starve_pH_5_sig = HtoA_snf5_Dex_v_Starve_pH_5[HtoA_snf5_Dex_v_Starve_pH_5$padj < 0.05,]
#filter log2FC over 1 (two-fold)
  WT_Dex_v_Starve_pH_7_sig = WT_Dex_v_Starve_pH_7[abs(WT_Dex_v_Starve_pH_7$log2FoldChange) > 1,]
  WT_DEX_v_starve_pH_5_sig = WT_Dex_v_Starve_pH_5[abs(WT_Dex_v_Starve_pH_5$log2FoldChange) > 1,]
  dQ_snf5_Dex_v_Starve_pH_5_sig = dQ_snf5_Dex_v_Starve_pH_5[abs(dQ_snf5_Dex_v_Starve_pH_5$log2FoldChange) > 1,]
  Null_Dex_v_Starve_sig = Null_Dex_v_Starve[abs(Null_Dex_v_Starve$log2FoldChange) > 1,]
  HtoA_snf5_Dex_v_Starve_pH_5_sig = HtoA_snf5_Dex_v_Starve_pH_5[abs(HtoA_snf5_Dex_v_Starve_pH_5$log2FoldChange) > 1,]
  
# Writing the list to CSV
write.table(as.data.frame(WT_Dex_v_Starve_pH_7_sig), file="WT_Dex_v_Starve_pH_7_sig.txt", sep ="\t")
write.table(as.data.frame(WT_DEX_v_starve_pH_5_sig), file="WT_DEX_v_starve_pH_5_sig.txt", sep ="\t")
write.table(as.data.frame(dQ_snf5_Dex_v_Starve_pH_5_sig), file="dQ_snf5_Dex_v_Starve_pH_5_sig.txt", sep ="\t")
write.table(as.data.frame(Null_Dex_v_Starve_sig), file="Null_Dex_v_Starve_sig.txt", sep ="\t")
write.table(as.data.frame(HtoA_snf5_Dex_v_Starve_pH_5_sig), file="HtoA_snf5_Dex_v_Starve_pH_5_sig.txt", sep ="\t")



#WT_Dex_v_Starve_pH_5$padj < 0.1, na.rm=TRUE)

