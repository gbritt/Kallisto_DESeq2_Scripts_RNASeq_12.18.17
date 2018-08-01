#DESeq2 LRT + heatmap

#when running this run 'prepare' by hand without running hte DESeq comman. Not sure if doing the LRT overwrites the Padj values from before or not. SO have a fresh start each time

library(pheatmap)
library(DESeq2)
library(RColorBrewer)
source('/Users/HoltLab/Documents/R/is.true.R') # function should be in directory

setwd("~/Documents/R/Scripts/Kallisto_DESeq2_Scripts_RNASeq_12.18.17/")

dds = readRDS('dds.rda') # for individual: dds_individual.rda
sampledata = readRDS('sampledata.rda') # for individual: sampledata_individual.rda
t2g = readRDS('t2g.rda')
DE_Genes_in_Starvation = readRDS("DE_Genes_in_Starvation_combined_master.rda")  # individual DE_Genes_in_Starvation_master_uncut.rda

results_path = '/Users/HoltLab/Documents/R/Experiments/RNA_seq_ADH2_Nacho_12.28.17/Results/DESeq2/DE_genes_combined/' # choose DE_genes_individual_comparisons if you want it

dQ_stv_v_dex = read.table(file = paste(results_path, "dQ_Starve_v_Dex.txt", sep =""), sep = "\t") #
Null_stv_v_dex = read.table(file = paste(results_path, "Null_Starve_v_Dex.txt", sep =""), sep = "\t") # 
WT_stv_v_dex = read.table(file = paste(results_path, "WT_Starve_v_Dex.txt", sep =""), sep = "\t") #
HtoA_stv_v_dex = read.table(file = paste(results_path, "HtoA_Starve_v_Dex.txt", sep =""), sep = "\t") #
WT_stv_v_dex_pH7 = read.table(file = paste(results_path, "WT_Starve_v_Dex_pH_7.txt", sep =""), sep = "\t") #


# dQ_WT_interaction_stv_v_dex = read.table(file = paste(results_path,"dQ_WT_interaction_stv_v_dex.txt", sep = ""), sep ="\t")
# HtoA_WT_interaction_stv_v_dex = read.table(file = paste(results,"HtoA_WT_interaction_stv_v_dex.txt", sep = ""), sep ="\t")
# Null_WT_interaction_stv_v_dex = read.table(file = paste(results,"Null_WT_interaction_stv_v_dex.txt", sep = ""), sep ="\t")
# dQ_HtoA_interaction_stv_v_dex = read.table(file = paste(results_path,"dQ_HtoA_interaction_stv_v_dex.txt", sep = ""), sep ="\t")


Cluster_Genes = c(rownames(WT_stv_v_dex),rownames(dQ_stv_v_dex),rownames(HtoA_stv_v_dex),rownames(Null_stv_v_dex),rownames(WT_stv_v_dex_pH7))
Cluster_Genes = c(rownames(WT_stv_v_dex),rownames(dQ_stv_v_dex),rownames(HtoA_stv_v_dex))

Cluster_Genes = unique(unlist(strsplit(Cluster_Genes, " "))) # all the DE genes

heatmap_genes = DE_Genes_in_Starvation[Cluster_Genes,] # list of all log2FC - should from from 'prepare' script  # for when you want to make a heatmap of log2FC
heatmap_genes = res_ordered[DEstarveGenes@rownames,] # gives you expression values             # For when you want to make a heatmap of expression values 
x <- scale(heatmap_genes) # scale and center columns
x <- t(scale(t(heatmap_genes))) # scale and center rows
x = na.omit(x)

y = x[,c(1:3)] # if you want to drop some samples plot with y

df = as.data.frame(c(rep("WT_Dex",3),rep("dQ_Dex",3),rep("HtoA_Dex",3),rep("Null_Dex",3),rep("WT_Stv_pH5",4),rep("dQ_Stv_pH5",4),rep("HtoA_Stv_pH5",4),rep("Null_Stv_pH5",4),rep("WT_Stv_pH7",3)))

df = as.data.frame(c(rep("WT_Dex",3),rep("dQ_Dex",3),rep("HtoA_Dex",3),rep("Null_Dex",3),rep("WT_Stv_pH5",3),rep("WT_Stv_pH7",3),rep("dQ_Stv_pH5",3),rep("HtoA_Stv_pH5",3),rep("Null_Stv_pH5",3))) # for when I use shortened input

col2 = as.data.frame(c(rep("Dex",12), rep("Starve",19 ))) # before I cut a bunch off
col2 = as.data.frame(c(rep("Dex",12),rep("Starve",15)))

df = cbind(df,col2)
colnames(df) = c("Condition", "Growth Condition")

rownames(df) <- colnames(heatmap_genes)


pl = pheatmap(x, cluster_rows=TRUE, show_rownames=F, 
         cluster_cols=F, show_colnames = T, border_color = "NA", fontsize_col = 20, , annotation_col = df, fontsize_row = 8), cutree_rows = 7) #or use heatmap genes
pheatmap(x, border_color = "grey60", annotation_col = df)

pl = pheatmap(x, border_color = "NA",show_rownames = F,cluster_cols = F,  annotation_col = df, cutree_rows = 4) # for gene expression


#y.order <- cbind(y[c(pl$tree_row[["order"]],pl$tree_col[["order"]]]),cluster=cutree(pl$tree_row, k=7)[pl$tree_row[["order"]]]) # when you clustered columns and rows # the c() might have been in the wrong place?

pl = pheatmap(x, border_color = "NA",show_rownames = F,cluster_cols = F,  annotation_col = df, cutree_rows = 14) # for gene expression

lbl <- cutree(pl$tree_row, 14) #need to change to same cutree_rows from pheatmap


x.order <- cbind(x[pl$tree_row[["order"]]],cluster=cutree(pl$tree_row, k=14)[pl$tree_row[["order"]]]) # when you clustered rows



#above gives genes log2FC and cluster
unique(x.order[,"cluster"]) #gives the cluster order
# so here it is 1 4 3 6 2 7 5

names_cluster1 = which(lbl==10) # find genes corresponding to first group, ...
names_cluster2 = which(lbl==6)
names_cluster3 = which(lbl==11)
names_cluster4 = which(lbl==12)
names_cluster5 = which(lbl==9)
names_cluster6 = which(lbl==13)
names_cluster7 = which(lbl==1)
names_cluster8 = which(lbl==8)
names_cluster9 = which(lbl==14)
names_cluster10 = which(lbl==7)
names_cluster11 = which(lbl==3)
names_cluster12 = which(lbl==2)
names_cluster13 = which(lbl==5)
names_cluster14 = which(lbl==4)

# 6 clusters for gene expression clustering!!!!!

names_cluster1 = which(lbl==3) # find genes corresponding to first group, ...
names_cluster2 = which(lbl==1)
names_cluster3 = which(lbl==4)
names_cluster4 = which(lbl==2)

cluster_list = list(names_cluster1,names_cluster2,names_cluster3,names_cluster4, names_cluster5, names_cluster6, names_cluster7, names_cluster8, names_cluster9, names_cluster10, names_cluster11, names_cluster12, names_cluster13, names_cluster14)
cluster_writer = function(cluster_list){ 
  setwd("/Users/HoltLab/Documents/R/Scripts/Kallisto_DESeq2_Scripts_RNASeq_12.18.17")
  library('gProfileR')
  library(knitr)
  library(kableExtra)
  library('magick')
  
  for(i in 1:length(cluster_list)){
    write.table(as.data.frame(cluster_list[i]), file= paste('Cluster_Genes_group',i,'.txt', sep = ""), sep ="\t")
    
    #Now do GO enrichment (from GO_Enrichment.R)
    gene_name = names(cluster_list[[i]])
    term.induced <- gprofiler(query=gene_name, organism="scerevisiae")
    term.induced <- term.induced[order(term.induced$p.value),]
    # term.induced$p.value
    results_path = '/Users/HoltLab/Documents/R/Experiments/RNA_seq_ADH2_Nacho_12.28.17/Results/DESeq2/'
    
    # Kable allows you to create pretty markdown/latex/PDF version of table - can worry about this later, for now it just prints well in R (won't work in function)
    kable(term.induced[1:10,c("term.name",
                              "term.size",
                              "query.size",
                              "overlap.size",
                              "recall",
                              "precision",
                              "p.value", 
                              "intersection")], 
    format.args=c(engeneer=TRUE, digits=3), caption="**Table: functional analysis wit gProfileR. ** ", keep_pdf = TRUE)
    
    
    write.table(term.induced, file= paste(results_path, 'Cluster_Genes_group',i, "_GO_enrichment", ".txt", sep = ""), sep ="\t")
  }
      
      
  
}


write.table(as.data.frame(names_cluster7), file= 'Genes_Group7_(5).txt', sep ="\t") # From here on out will name based on how this cluster falls in heatmap rather than how the clustering chose them

  
