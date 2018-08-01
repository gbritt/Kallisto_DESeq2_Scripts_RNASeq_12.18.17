
df = as.data.frame(WT_Starve_v_Dex_Volcano)



Volcano1 = function(sample_name){
  library('ggrepel')
  library('ggplot2')
  df = as.data.frame(sample_name)
  genes = df # I am lazy and leaving this
  genes$Significant <- ifelse(df$padj < 0.05 & abs(genes$log2FoldChange) > 1, "padj < 0.05", "Not Sig")
  genes$Gene = genes$ext_gene
  
  ggplot(genes, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = Significant)) +
    ylim(0,60) + #setting limits so all volcanos have the same axes
    xlim(-8,8) +
    scale_color_manual(values = c( "grey", "red")) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") +
    geom_text_repel(
      data = subset(genes, padj < 0.0000001 & abs(log2FoldChange) > 3),
      aes(label = Gene),
      size = 3,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")
    )
  
  
}



# Orignal Try - can use if you want to cut off the top ----
df = as.data.frame(test)
df$shape <- ifelse(df$padj < 0.0000000000003, "triangle", "circle")
df$padj[df$padj<0.0000000000003] = 0.0000000000003

df$shape[(abs(df$log2FoldChange) > 3)] <- "triangle"
df$log2FoldChange[df$log2FoldChange >  3] <- 3
df$log2FoldChange[df$log2FoldChange < -3] <- 3

sp = ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj), shape=shape), label = rownames(df)) +
  geom_point(alpha=0.5, size=4) +
  theme(legend.position = "none") +
  xlim(c(-1.5, 1.5)) + ylim(c(0, -log10(0.0000000000003))) +
  xlab("log2 fold change") + ylab("-log10 p-value")

sp = ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj)), label = rownames(df)) +
  geom_point(alpha=0.5, size=4) +
  theme(legend.position = "none") +
  xlim(c(-1.5, 1.5))  +
  xlab("log2 fold change") + ylab("-log10 p-value")