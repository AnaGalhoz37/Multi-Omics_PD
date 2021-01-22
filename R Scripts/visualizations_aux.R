# Volcano plot
volcano_dist_plot <- function(res_matrix,alpha,title){
  col_data <- densCols(res_matrix$'log2FoldChange', 
                       -log10(res_matrix$'pvalue'))
  graphics::plot(res_matrix$'log2FoldChange',-log10(res_matrix$'padj'),col = col_data,
                 main = title,
                 xlab = "Effect size: log2(fold-change)",ylab = "-log10(adjusted p-value)",
                 pch=20,cex=0.5,bty = "n",xlim = c(-1,1),ylim = c(0,3))
  # with(subset(data.frame(res_matrix), padj<0.5), 
  #      points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
  with(subset(data.frame(res_matrix), padj>0.1), 
       points(log2FoldChange, -log10(padj), pch=20, col=gray(level = 0.8,alpha = 0.7), cex=0.6))
  with(subset(data.frame(res_matrix), padj<0.1 & log2FoldChange > 0), 
       points(log2FoldChange, -log10(padj), pch=20, col="#ff6961", cex=0.5))
  with(subset(data.frame(res_matrix), padj<0.1 & log2FoldChange < 0), 
       points(log2FoldChange, -log10(padj), pch=20, col="#6A9EFC", cex=0.5))
  # with(subset(data.frame(res_matrix), padj<alpha),
  #      points(log2FoldChange, -log10(padj), pch=20, col = "green", cex=0.5))
  abline(v=0)
  #abline(v=c(-1,1), col="brown")
  abline(h=-log10(0.1), col="brown")
}

# Pairs between miRNA and Gene symbols as an example 
# Input: vectors of Log2fc of miRNA and Gene Symbol in correct order 
paired_miRNA_gene <- function(miRNA, Gene_Symbol,title){
  x = data.frame(miRNA, Gene_Symbol)
  xx<-stack(x) # restructures data for stripchart function
  par(las=1) # horizontal axis labels
  with(xx,  stripchart(values~ind, xlim=c(0,3),pch=19, 
                       main=title,ylab="Pairs miRNA-Gene", vertical=TRUE, 
                       col=c("red", "green"),axes = FALSE))
  axis(1)
  axis(2)
  abline(h = 0, col = "grey",lty = 2)
  apply(x,1,lines, col="black") 
}

# Boxplots 
boxplot_validation <- function(data_validate,col_data,omics_validate){
  for (i in 1:length(omics_validate)) {
    omic <- paste0('^',omics_validate[i],'$')
    data.omic <- data.frame(t(data_validate[grep(omic,rownames(data_validate)),]))
    if(length(colnames(data.omic))!=0){
      data.omic$'ID' <- rownames(data.omic)
      data.omic$'condition' <- col_data[,1]
      data.omic$'batch' <- col_data[,2]
      names(data.omic) <- c("omic","ID","condition","batch")
      temp_plot = ggplot(data.omic,aes(condition,omic)) +
        geom_boxplot(aes(fill = factor(condition)),alpha = 0.7,outlier.shape = NA) +
        geom_point(aes(color = factor(batch)),
                   position = position_jitter(width = 0.1, height = 0))+
        #geom_text(aes(label = ID)) +
        labs(x = "",
             y="Normalized Counts",
             title = omics_validate[i]) +
        #geom_jitter(position=position_jitter(width=.2, height=0)) +
        scale_fill_manual(breaks = c("treated","untreated"),
                          values = c("#b19cd9","#ffb347")) +
        scale_color_manual(breaks = c("1","2"),
                           values=c("red", "black"),
                           name = "Batch",
                           labels = c("Validation","Discovery")) +
        scale_x_discrete(labels = c("PD","CTR")) +
        theme_classic() +
        theme(axis.text = element_text(size=12),
              title = element_text(size=12)) +
        guides(fill = FALSE)
      ggsave(temp_plot,path = "Plots_from_code",file = paste0("box-plot_",omics_validate[i],".pdf"),
             width=4, height = 4.5, units="in")
    }
  }
}

