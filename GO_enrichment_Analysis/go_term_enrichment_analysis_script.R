
library(readr)
library(dplyr)
library(tidyr)
library(magrittr)
#args[1] gives the name of the module we are analyzing

moduleGOTermCounts <- function(module)
{
  fin <- read_tsv(paste0("E:/AshtonData/MPI_all_data/Final_results_hopefully/", module, "/submodule_go_terms.tsv"))
  count_d <- fin %>% count(GO_ID) 
  count_d %>% mutate("p-vals" = pValReporting(count_d, module))
}

#Read in what you want to compare it to

pValReporting <- function(counts_list, module_name)
{
  file_name <- paste0("E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/", module_name, ".txt")
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

go_eval <- function(z){
count_f <- moduleGOTermCounts(z)

write_tsv(count_f, paste0("E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/", z, "_go-p-vals.tsv"))
}

go_eval("NIATv7_g15400")
go_eval("NIATv7_g40277")
go_eval("NIATv7_g04345")
go_eval("NIATv7_g41243")
go_eval("NIATv7_g19262")
go_eval("NIATv7_g41088")
go_eval("NIATv7_g30725")
go_eval("NIATv7_g27755")
go_eval("NIATv7_g03410")
go_eval("NIATv7_g05680")
go_eval("NIATv7_g19992")
go_eval("NIATv7_g26487")
go_eval("NIATv7_g33088")
go_eval("NIATv7_g30738")
go_eval("NIATv7_g32572")
go_eval("NIATv7_g33798")
go_eval("NIATv7_g07696")
go_eval("NIATv7_g12711")
go_eval("NIATv7_g16189")
go_eval("NIATv7_g17075")
go_eval("NIATv7_g18001")
go_eval("NIATv7_g19088")
go_eval("NIATv7_g21131")
go_eval("NIATv7_g28241")
go_eval("NIATv7_g29978")
go_eval("NIATv7_g33798")
go_eval("NIATv7_g34810")