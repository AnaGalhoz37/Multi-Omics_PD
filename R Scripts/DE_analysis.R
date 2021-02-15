# Define count data for PD and CTR patients
count.PD <- demo.count.file[,1:10]
count.CTR <- demo.count.file[,11:20]

# Column data with the condition corresponding to each ID in count matrix
col.count.data <- col_data(demo.count.file,count.PD,count.CTR)

# -> Creation of matrix with the count information 
mRNA.count.data <- matrix_raw_data(demo.count.file,col.count.data)

# -> Creation of DESeq2 Data Matrix
dds <- DESeqDataSetFromMatrix(countData = mRNA.count.data,
                              colData = col.count.data,
                              design = ~ condition)

# Eliminate the zero counts from mRNA data (pre-filter)
mRNA.count.tokeep <- rowSums(counts(dds)) > 0
mRNA.count.data.new <- dds[mRNA.count.tokeep,]

# Filter using TCGA package (based on quantile info)
dds <- estimateSizeFactors(dds)
mRNA.count.data.new.temp <- TCGAanalyze_Filtering(assay(dds), method = "quantile")
mRNA.count.data.new <- DESeqDataSetFromMatrix(countData = mRNA.count.data.new.temp,
                                              colData = col.count.data,
                                              design = ~ condition)

# Check which genes were removed in the filtering stage
house.keeping.genes <- assay(dds)[!rownames(dds) %in% rownames(assay(mRNA.count.data.new)),]
house.keeping.genes.map <- mRNA_map(house.keeping.genes,HGNC_genes,Ensembl_genes)

# Factor levels
mRNA.count.data.new$condition <- factor(mRNA.count.data.new$condition,
                                        levels = c("untreated","treated"))


# DE analysis
DE.mRNA.new <- DESeq(mRNA.count.data.new)
# default: p-value 0.1
res.DE.mRNA.new <- results(DE.mRNA.new,
                           name = "condition_treated_vs_untreated")

# Shrinkage of LFC for default p-value
res.DE.mRNA.new_shrink_apeglm <- lfcShrink(DE.mRNA.new,
                                   coef = "condition_treated_vs_untreated",
                                   type = "apeglm")
res.DE.mRNA.new_shrink_ashr <- lfcShrink(DE.mRNA.new,
                                   coef = "condition_treated_vs_untreated",
                                   type = "ashr")
# Ordering results from DESeq2 by p-value with default p-value
res.DE.mRNA.new.ordered <- res.DE.mRNA.new[order(res.DE.mRNA.new$pvalue),] 

# Get normalized data without any transformation (VSD or rlog)
normed_counts_mRNA <- counts(DE.mRNA.new,normalized = TRUE)

# Join normed counts with DESeq2 results
mRNA_counts_stats <- cbind(res.DE.mRNA.new,normed_counts_mRNA)

###########
# DE analysis 
#########
# Up and down regulated mRNA 
# -> not removed control
mRNA_up_reg <- mRNA_counts_stats[mRNA_counts_stats$log2FoldChange>1,]
mRNA_down_reg <- mRNA_counts_stats[mRNA_counts_stats$log2FoldChange<(-1),]

# Map up and down regulated mRNA with database (Ensembl - Gene name)
# -> not removed
mRNA_up_mapped <- mRNA_map(mRNA_up_reg,HGNC_genes,Ensembl_genes)
mRNA_down_mapped <- mRNA_map(mRNA_down_reg,HGNC_genes,Ensembl_genes)

# -> mapped
mRNA_up_Map <- gene_name(mRNA_up_mapped)
gene_symbol_up <- length(unique(mRNA_up_Map$`Gene Name`))
mRNA_down_Map <- gene_name(mRNA_down_mapped)
gene_symbol_down <- length(unique(mRNA_down_Map$`Gene Name`))

# mRNA with Gene Symbol
mRNA_Map <- data.frame(rbind(mRNA_up_Map,mRNA_down_Map))

# statistical genes
mRNA_Map_stat <- mRNA_Map[mRNA_Map$padj < 0.1 & !is.na(mRNA_Map$padj),]
