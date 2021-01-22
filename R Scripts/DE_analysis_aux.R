# auxiliary functions for the DEx analysis

# Create column data for diseased and healthy patients
col_data <- function(data_all,data_PD,data_CTRL){
  data_out <- as.data.frame(cbind(c("",colnames(data_all)),
                                  c("condition",rep("treated",ncol(data_PD)),
                                    rep("untreated",ncol(data_CTRL)))))
  colnames(data_out) <- as.character(unlist(data_out[1,]))
  data_out <- data_out[-1,]
  rownames(data_out) <- data_out[,1]
  data_out[,1] <- NULL
  return(data_out)
}

# Create matrix data for diseased and healthy patients
matrix_raw_data <- function(data_set,col_data){
  data_out <- data_set
  data_out <- data_out[,rownames(col_data)]
  return(data_out)
}

mRNA_HGNC <- function(data_set,data_HGNC){
  data_join <- merge(data_set,data_HGNC[,c("Ensembl gene ID","Approved symbol")],
                     by.x = 0,
                     by.y = "Ensembl gene ID",all.x=TRUE)
  nr_ids <- ncol(data_join)
  colnames(data_join)[nr_ids] <- "Gene Name" 
  colnames(data_join)[1] <- "Ensembl gene ID"
  return(data_join)
}

mRNA_Ensembl <- function(data_set,data_Ensembl){
  data_join <- merge(data_set,data_Ensembl,
                     by.x = 0,
                     by.y = "Gene stable ID",all.x=TRUE)
  nr_ids <- ncol(data_join)
  colnames(data_join)[nr_ids] <- "Gene Name" 
  colnames(data_join)[1] <- "Ensembl gene ID"
  return(data_join)
}

mRNA_map <- function(data_set,data_HGNC,data_Ensembl){
  mRNA_data_HGNC <- mRNA_HGNC(data_set,data_HGNC)
  mRNA_data_Ensembl <- mRNA_Ensembl(data_set,data_Ensembl)
  mRNA_data_HGNC_noName <- NA_gene_name(mRNA_data_HGNC)
  mRNA_data_HGNC_Name <- With_gene_name(mRNA_data_HGNC,mRNA_data_HGNC_noName)
  data_aux <- mRNA_data_Ensembl[mRNA_data_Ensembl$"Ensembl gene ID" 
                                %in% mRNA_data_HGNC_noName$"Ensembl gene ID",]
  data_aux <- rbind(mRNA_data_HGNC_Name,data_aux)
  return(data_aux)
}

# Identification of Ensembl gene ID's with no Gene Symbol
NA_gene_name <- function(data_PD_CTRL){
  data_NA <- data_PD_CTRL[is.na(data_PD_CTRL$"Gene Name"),]
}

With_gene_name <- function(data_set_all,data_na){
  data_with <- setdiff(data_set_all$"Ensembl gene ID", data_na$"Ensembl gene ID")
  data_out <- data_set_all[data_set_all$"Ensembl gene ID" %in% data_with,]
}

gene_name <- function(data_PD_CTRL){
  data_NA <- data_PD_CTRL[!is.na(data_PD_CTRL$"Gene Name"),]
  return(data_NA)
}
