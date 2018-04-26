#Generic version of the script to make the p-values.
library(readr)
library(dplyr)
library(tidyr)
library(magrittr)
#args[1] gives the name of the module we are analyzing

moduleGOTermCounts <- function(path_to_gos, path_to_counts)
{
  fin <- read_tsv(path_to_gos)
  count_d <- fin %>% count(GO_ID) 
  count_d %>% mutate("p-vals" = pValReporting(count_d, path_to_counts))
}

#Read in what you want to compare it to

pValReporting <- function(counts_list, path_to_perms)
{
  file_name <- path_to_perms
  f_in <- read_tsv(file_name)
  apply(counts_list, 1, function(x) pValCalc(x, f_in))
}

pValCalc <- function(counts_list, permuted_data)
{
  count <- as.numeric(counts_list[2])
  go_id <- counts_list[1]
  permuted_data %>% filter(GO_ID == go_id) -> sto
  unlist(as.numeric(sto[-1]),use.names = F) -> sto
  #Lookup the value
  val <- sum(sto >= count)
  val/1000
}

go_eval <- function(path_to_gos, path_to_counts, out_path){
  count_f <- moduleGOTermCounts(path_to_gos,path_to_counts )
  
  write_tsv(count_f, out_path)
}
setwd("E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/")
go_eval("MYC2_MA_submodule/myc2_go_terms.tsv","MYC2_MA_submodule/MYC2_p_counts.txt", "MYC2_MA_submodule/myc_2_go_p-vals.tsv")


go_eval("MYC2_MA_submodule/myc2_go_terms.tsv","MYC2_MA_submodule/MYC2_p_counts.txt", "MYC2_MA_submodule/myc_2_go_p-vals.tsv")

setwd("E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/MYC2_RNA_submodule/")
go_eval("MYC2_RNA_go_terms.tsv", "MYC2_RNA_submodule_p_counts.txt", "myc_2_go_p-vals.tsv")


setwd("E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/WRKY3_MA/")
go_eval("go_list.txt", "WRKY3_MA_submodule_p_counts.txt", "wrky3_go_p-vals.tsv")


