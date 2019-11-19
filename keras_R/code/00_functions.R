library(seqinr)
library(stringr)
library(reshape2)
library(data.table)
library(dplyr)
library(gtools)
library(Biostrings)
options(datatable.fread.input.cmd.message=FALSE)

# Lots of misc

importForDeep <- function(sample, path){
 
  editor <- "ABE"
  saveForEval <- c("chr21", "chr22", "chrX")
  
  # Import sequence fasta files
  dir_seq <- paste0(path,"/fastas/", editor, "-sequence")
  
  seq_fastas <- list.files(dir_seq, pattern = paste0(sample), full.names = TRUE)
  seq_fastas <- seq_fastas[grepl(".gz", seq_fastas) & !grepl(paste(saveForEval,collapse="|"), seq_fastas)]
  
  # Import data into a list
  fasta_input <- sapply(seq_fastas, function(faf) unlist(read.fasta(faf, as.string = TRUE)))
  
  # Get all sequence attributes
  seq_names <- sapply(fasta_input, names) %>% unlist() %>% unname()
  meta <- data.frame(stringr::str_split_fixed(seq_names, pattern = "_", n = 7), stringsAsFactors = FALSE)
  colnames(meta) <- c("Library", "chr_pos", "strand", "editRate", "coverage",  "annotation", "gene")
  meta$chr <- stringr::str_split_fixed(meta$chr_pos, "-", 2)[,1]
  seqs <- unlist(fasta_input) %>% unname()
  
  # Process meta to convert characters to numeric
  meta$editRate <- as.numeric(meta$editRate); meta$coverage <- as.numeric(meta$coverage)
  rownames(meta) <- NULL
  meta$sequence <- seqs
  meta$isEdited <- meta$editRate > 0.05
  boo <- meta$editRate > 0.05 | meta$editRate == 0
  return(meta[boo,])
}

importForDeep_wStructure <- function(sample, editor){
  
  saveForDeepLift <- c("chr21", "chr22", "chrX")
  
  # Import sequence fasta files
  dir_seq <- paste0("../../CRISPR-",editor,"-RNA-DATA/fastas/", editor, "-sequence")
  dir_struct <- paste0("../../CRISPR-",editor,"-RNA-DATA/fastas/", editor, "-secondary")
  
  seq_fastas <- list.files(dir_seq, pattern = paste0(sample), full.names = TRUE)
  seq_fastas <- seq_fastas[grepl(".gz", seq_fastas) & !grepl(paste(saveForDeepLift,collapse="|"), seq_fastas)]
  struct_files <- list.files(dir_struct, pattern = paste0(sample), full.names = TRUE)
  struct_files <- struct_files[grepl(".gz", struct_files) & !grepl(paste(saveForDeepLift,collapse="|"), struct_files)]
  
  # Import data into a list
  fasta_input <- sapply(seq_fastas, function(faf) unlist(read.fasta(faf, as.string = TRUE)))
  
  # Process structure data
  struct_input <- sapply(struct_files, function(faf) fread(cmd = paste0("zcat < ", faf), header = FALSE)[[1]][c(FALSE, TRUE)]) %>% unlist() %>% unname()
  
  # Get all sequence attributes
  seq_names <- sapply(fasta_input, names) %>% unlist() %>% unname()
  meta <- data.frame(stringr::str_split_fixed(seq_names, pattern = "_", n = 7), stringsAsFactors = FALSE)
  colnames(meta) <- c("Library", "chr_pos", "strand", "editRate", "coverage",  "annotation", "gene")
  meta$chr <- stringr::str_split_fixed(meta$chr_pos, "-", 2)[,1]
  seqs <- unlist(fasta_input) %>% unname()
  
  # Process meta to convert characters to numeric
  meta$editRate <- as.numeric(meta$editRate); meta$coverage <- as.numeric(meta$coverage)
  rownames(meta) <- NULL
  meta$sequence <- seqs
  meta$structure <- struct_input
  meta$isEdited <- meta$editRate > 0.05
  boo <- meta$editRate > 0.05 | meta$editRate == 0
  return(meta[boo,])
}

make_one_hot_eight_channel <- function(x) {
  
  # Split into each individual matrix
  seq_mat <-str_split_fixed( x[,1],  "", 101)
  struct_mat <- str_split_fixed(x[,2],  "", 101)
  
  # Encode a fifth channel for sequence -- editing event
  seq_mat[,51] <- "Z"
  
  # Rehsape into nseqs x sequence length x nucleotide
  rm <- rbind(reshape2::melt(seq_mat), reshape2::melt(struct_mat)); rm$new <- 1
  aa <- acast(rm, Var1  ~  Var2  ~  value, fill = 0, value.var = "new")  
  aa
}


make_one_hot_five_channel_string <- function(x) {
  
  # Split sequences into individual characters
  spm <- str_split_fixed(x, "", 101)
  
  # Encode a fifth channel
  spm[,51] <- "Z"
  
  # Rehsape into nseqs x sequence length x nucleotide
  rm <- reshape2::melt(spm); rm$new <- 1
  aa <- acast(rm, Var1  ~  Var2  ~  value, fill = 0, value.var = "new")  
  aa
}
