### RNA decomposition

# 2nd batch: with best (cpm)
data_2 <- mRNA_map(mRNA.count.data.new,HGNC_genes,Ensembl_genes)
data_2 <- na.omit(data_2)
rownames(data_2) <- make.unique(data_2$`Gene Name`)
data_2 <- data_2[,-c(1,ncol(data_2))]
data_2_tpm <- Counts_to_tpm(data_2,rep(50,nrow(data_2)))
data_2_PD_tpm <- data_2_tpm[,names(data_2) %in% colnames(count.PD)]
data_2_CTRL_tpm <- data_2_tpm[,names(data_2) %in% colnames(count.CTR)]

# -> xcell
deconv_xcell_PD_2 <- immunedeconv::deconvolute(data_2_PD_tpm,"xcell")
deconv_xcell_PD_2$mean.PD <- apply(deconv_xcell_PD_2[,-1],1,mean)
deconv_xcell_CTRL_2 <- immunedeconv::deconvolute(data_2_CTRL_tpm,"xcell")
deconv_xcell_CTRL_2$mean.CTRL <- apply(deconv_xcell_CTRL_2[,-1],1,mean)
deconv_xcell_2 <- immunedeconv::deconvolute(data_2_tpm,"xcell")
deconv_xcell_2$log2FC <- log2(deconv_xcell_PD_2$mean.PD/deconv_xcell_CTRL_2$mean.CTRL)
deconv_xcell_2$p.value <- rep(1,nrow(deconv_xcell_2))
for(i in 1:nrow(deconv_xcell_2)){
  deconv_xcell_2[i,]$p.value <- t.test(deconv_xcell_PD_2[i,-1],deconv_xcell_CTRL_2[i,-1])$p.value
}
deconv_xcell_2 <- deconv_xcell_2[order(deconv_xcell_2$p.value,decreasing = T),]
deconv_xcell_2$p.adjusted <- p.adjust(deconv_xcell_2$p.value,method = "BH")

# -> mcp_counter
deconv_mcpcounter_PD_2 <- immunedeconv::deconvolute(data_2_PD_tpm,"mcp_counter")
deconv_mcpcounter_PD_2$mean.PD <- apply(deconv_mcpcounter_PD_2[,-1],1,mean)
deconv_mcpcounter_CTRL_2 <- immunedeconv::deconvolute(data_2_CTRL_tpm,"mcp_counter")
deconv_mcpcounter_CTRL_2$mean.CTRL <- apply(deconv_mcpcounter_CTRL_2[,-1],1,mean)
deconv_mcpcounter_2 <- immunedeconv::deconvolute(data_2_tpm,"mcp_counter")
deconv_mcpcounter_2$log2fc <- log2(deconv_mcpcounter_PD_2$mean.PD/deconv_mcpcounter_CTRL_2$mean.CTRL)
deconv_mcpcounter_2$p.value <- rep(1,nrow(deconv_mcpcounter_2))
for(i in 1:nrow(deconv_mcpcounter_2)){
  deconv_mcpcounter_2[i,]$p.value <- t.test(deconv_mcpcounter_PD_2[i,-1],deconv_mcpcounter_CTRL_2[i,-1])$p.value
}
deconv_mcpcounter_2 <- deconv_mcpcounter_2[order(deconv_mcpcounter_2$p.value,decreasing = T),]
deconv_mcpcounter_2$p.adjusted <- p.adjust(deconv_mcpcounter_2$p.value,method = "BH")

# --------------- auxiliar function
Counts_to_tpm <- function(counts, featureLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- featureLength
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen)
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}
