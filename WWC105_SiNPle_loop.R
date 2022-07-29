library(data.table)
library(purrr)
library(dplyr)
library(qpcR)
library(UpSetR)
library(ggplot2)
library(janitor)
library(openxlsx)

file_list <- Sys.glob("*variants*.txt")
head(file_list)

for (i in seq_along(file_list)) 
  {
  filename <- file_list[[i]]

  # read a file 
  sinple_variant <- read.table(filename, sep="\t", header=FALSE, col.names = paste("V", 1:50),fill=TRUE)
  head(sinple_variant)
  #sinple_variant <- sinple_variant[,which(is.na(sinple_variant[1,]))]
  #head(sinple_variant)

  alpha <- fread("alpha.txt", header=TRUE)
  #head(alpha)
  names(alpha)[1] <- "V.2"
  #head(alpha)

  common_mutation <- sinple_variant %>% filter(V.2 %in% alpha$V.2)

  name_list <- unlist(strsplit(filename, "_",2))
  
  file1 <- paste("common", "alpha" ,name_list[1], name_list[2], sep="_")
  file2 <- paste(file1, "xlsx", sep=".")

  write.xlsx(common_mutation, file2)
  }