# load necessary data

# R files
load(demo.count.file.RData)
load(MTI_mirTarBase.RData)
load(Ensembl_genes.RData)

dir_use <- "Data/"

# Get names of file lists
get_res <- function(read_file,file,min,mx){
  res_file <- lapply(read_file, summary)
  names(res_file) <- substr(file,min,mx)
  return(res_file)
}

# txt files
txt_files <- dir(dir_use, pattern = ".txt", full.names = TRUE)
read_txt <- lapply(txt_files,read.delim)

# xls files
excel_files <- dir(dir_use, pattern = ".xls",full.names = TRUE)
read_excel <- lapply(excel_files, read_excel)

# csv files
csv_files <- dir(dir_use, pattern = ".csv",full.names = TRUE)
read_csv <- lapply(csv_files, read_csv)
