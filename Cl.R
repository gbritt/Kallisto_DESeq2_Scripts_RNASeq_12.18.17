# cluster profiler
library('gProfileR')
library(knitr)
library(kableExtra)
library('magick')


a = names(names1)

GO_enrichment(name){
  term.induced <- gprofiler(query=name, organism="scerevisiae")
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
  
  write.table(term.induced, file= paste(results_path,"GO_enrichment" ,".txt", sep = ""), sep ="\t")
  
}


res.DESeq2.df <- na.omit(data.frame(res.DESeq2))
induced.sign <- rownames(res.DESeq2.df)[res.DESeq2.df$log2FoldChange >= 2 &  res.DESeq2.df$padj < alpha]
# head(induced.sign)
# names(term.induced)


