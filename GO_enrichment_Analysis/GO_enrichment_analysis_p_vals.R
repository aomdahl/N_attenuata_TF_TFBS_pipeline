
library(readr)
library(dplyr)
library(tidyr)
library(magrittr)
rna_module_sizes <- data.frame("Module" =c("turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink", "magenta", "purple", "greenyellow", "tan", "salmon", "cyan", "midnightblue", "lightcyan", "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise", "darkgrey", "orange", "darkorange", "white", "skyblue", "saddlebrown", "paleturquoise", "steelblue", "violet", "gre60"), "Size" = c(999, 982, 908, 846, 784, 749, 667, 604, 578, 526, 473, 279, 273, 148, 131, 119, 101, 74, 73, 72, 63, 58, 54, 51, 47, 45, 40, 38, 35, 35, 34, 114))
MA_module_sizes <- data.frame("Module" = c("turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink", "magenta", "purple", "greenyellow", "tan", "salmon", "cyan", "midnightblue", "lightcyan", "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise", "darkgrey", "orange", "darkorange", "white", "skyblue", "saddlebrown", "steelblue", "paleturquoise", "violet", "darkolivegreen", "darkmagenta", "grey60"), "Size" = c(946, 829, 790, 780, 677, 577, 552, 432, 379, 370, 362, 337, 328, 245, 213, 206, 173, 156, 128, 127, 119, 113, 112, 110, 106, 106, 100, 97, 77, 70, 66, 58, 53, 206))
setwd("E:/AshtonData/MPI_all_data/GO_enrichment_analysis")
ma_go_terms <- read_tsv("ma_go_genes_modules.tsv")
ma_go_terms %>% count(Module, GO_ID) %>% arrange(desc(n)) -> MA_go_term_counts

setwd("E:/AshtonData/MPI_all_data/GO_enrichment_analysis")
RNA_go_terms <- read_tsv("rna_go_genes_modules.tsv")
RNA_go_terms %>% count(Module, GO_ID) %>% arrange(desc(n)) -> RNA_go_term_counts

inner_join(RNA_go_term_counts, rna_module_sizes, by = c("Module")) %>% rename("ModuleSize" = Size) -> rna_enrichment_sizes
inner_join(MA_go_term_counts, MA_module_sizes, by = c("Module")) %>% rename("ModuleSize" = Size) -> MA_enrichment_sizes
write_tsv(MA_enrichment_sizes, "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_enrichment_sizes.tsv")

#Read in what you want to compare it to

pValReportingRNA <- function(module_name, go_id)
{
  file_name <- paste0("E:/AshtonData/MPI_all_data/GO_enrichment_analysis/RNA_permutations/rna_", module_name, ".txt")
  f_in <- read_tsv(file_name)
  f_in %>% filter(GO_ID == go_id) -> sto
  unlist(as.numeric(sto[-1]),use.names = F) -> sto
  #Lookup the value
  lookup_val <- (RNA_go_term_counts %>% filter(Module == module_name, GO_ID == go_id))$n
  val <- sum(sto >= lookup_val)
  val/1000
}
#only doing the ones with over 50 enriched genes.
c() -> p_vals
c(p_vals, pValReportingRNA("turquoise", "GO:0016021")) -> p_vals
c(p_vals, pValReportingRNA("blue", "GO:0016021")) -> p_vals
c(p_vals, pValReportingRNA("red", "GO:0016021")) -> p_vals
c(p_vals, pValReportingRNA("yellow", "GO:0003735")) -> p_vals
c(p_vals, pValReportingRNA("black", "GO:0016021")) -> p_vals
c(p_vals, pValReportingRNA("brown", "GO:0016021")) -> p_vals
c(p_vals, pValReportingRNA("green", "GO:0016021")) -> p_vals
c(p_vals, pValReportingRNA("magenta", "GO:0016021")) -> p_vals
c(p_vals, pValReportingRNA("yellow", "GO:0006412")) -> p_vals
c(p_vals, pValReportingRNA("yellow", "GO:0016021")) -> p_vals
c(p_vals, pValReportingRNA("pink", "GO:0016021")) -> p_vals
c(p_vals, pValReportingRNA("purple", "GO:0016021")) -> p_vals
c(p_vals, pValReportingRNA("black", "GO:0009941")) -> p_vals
c(p_vals, pValReportingRNA("black", "GO:0009535")) -> p_vals
c(p_vals, pValReportingRNA("green", "GO:0055114")) -> p_vals
c(p_vals, pValReportingRNA("greenyellow", "GO:0016021")) -> p_vals
c(p_vals, pValReportingRNA("blue", "GO:0005524")) -> p_vals
c(p_vals, pValReportingRNA("turquoise", "GO:0005524")) -> p_vals
c(p_vals, pValReportingRNA("yellow", "GO:0022625")) -> p_vals
c(p_vals, pValReportingRNA("red", "GO:0005524")) -> p_vals
c(p_vals, pValReportingRNA("black", "GO:0009570")) -> p_vals
c(p_vals, pValReportingRNA("yellow", "GO:0005524")) -> p_vals
c(p_vals, pValReportingRNA("black", "GO:0019288")) -> p_vals
c(p_vals, pValReportingRNA("black", "GO:0006098")) -> p_vals
c(p_vals, pValReportingRNA("black", "GO:0006364")) -> p_vals
c(p_vals, pValReportingRNA("magenta", "GO:0005524")) -> p_vals
c(p_vals, pValReportingRNA("tan", "GO:0016021")) -> p_vals
c(p_vals, pValReportingRNA("turquoise", "GO:0005886")) -> p_vals
c(p_vals, pValReportingRNA("green", "GO:0046872")) -> p_vals
c(p_vals, pValReportingRNA("greenyellow", "GO:0005524")) -> p_vals
c(p_vals, pValReportingRNA("black", "GO:0010027")) -> p_vals
c(p_vals, pValReportingRNA("black", "GO:0009507")) -> p_vals
c(p_vals, pValReportingRNA("salmon", "GO:0016021")) -> p_vals
c(p_vals, pValReportingRNA("black", "GO:0055114")) -> p_vals
c(p_vals, pValReportingRNA("green", "GO:0020037")) -> p_vals
c(p_vals, pValReportingRNA("black", "GO:0010207")) -> p_vals
c(p_vals, pValReportingRNA("brown", "GO:0055114")) -> p_vals
c(p_vals, pValReportingRNA("turquoise", "GO:0055114")) -> p_vals
c(p_vals, pValReportingRNA("red", "GO:0005886")) -> p_vals
c(p_vals, pValReportingRNA("red", "GO:0055114")) -> p_vals
c(p_vals, pValReportingRNA("turquoise", "GO:0005634")) -> p_vals
c(p_vals, pValReportingRNA("yellow", "GO:0005634")) -> p_vals
c(p_vals, pValReportingRNA("blue", "GO:0008152")) -> p_vals
c(p_vals, pValReportingRNA("purple", "GO:0005524")) -> p_vals
c(p_vals, pValReportingRNA("green", "GO:0005576")) -> p_vals
c(p_vals, pValReportingRNA("turquoise", "GO:0005829")) -> p_vals
c(p_vals, pValReportingRNA("black", "GO:0019252")) -> p_vals
c(p_vals, pValReportingRNA("greenyellow", "GO:0005634")) -> p_vals
c(p_vals, pValReportingRNA("yellow", "GO:0022627")) -> p_vals


rna_enrichment_sizes %>% filter(n >= 50) %>% mutate("p-vals" = p_vals) -> rna_enrichment_with_pvals


pValReportingMA <- function(module_name, goid)
{
  file_name <- paste0("E:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_permutations/MA_", module_name, "_.txt")
  f_in <- read_tsv(file_name)
  f_in %>% filter(GO_ID == goid) -> sto
  unlist(as.numeric(sto[-1]),use.names = F) -> sto
  #Lookup the value
  lookup_val <- (MA_go_term_counts %>% filter(Module == module_name, GO_ID == goid))$n
  val <- sum(sto >= lookup_val)
  val/1000
}
MA_p_vals <- c()

c(MA_p_vals, pValReportingMA("turquoise", "GO:0016021")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("blue", "GO:0016021")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("green", "GO:0016021")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("yellow", "GO:0016021")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("brown", "GO:0003735")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("brown", "GO:0006412")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("red", "GO:0016021")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("black", "GO:0016021")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("brown", "GO:0016021")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("tan", "GO:0016021")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("magenta", "GO:0016021")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("red", "GO:0009941")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("pink", "GO:0016021")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("red", "GO:0009535")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("greenyellow", "GO:0016021")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("red", "GO:0009570")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("yellow", "GO:0055114")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("brown", "GO:0022625")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("turquoise", "GO:0005524")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("turquoise", "GO:0055114")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("purple", "GO:0016021")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("grey60", "GO:0016021")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("salmon", "GO:0016021")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("red", "GO:0006364")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("red", "GO:0019288")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("blue", "GO:0055114")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("red", "GO:0010027")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("red", "GO:0006098")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("turquoise", "GO:0008152")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("black", "GO:0005524")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("blue", "GO:0005524")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("brown", "GOvim ma_:0005524")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("brown", "GO:0005634")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("red", "GO:0009507")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("red", "GO:0055114")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("green", "GO:0055114")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("red", "GO:0010207")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("green", "GO:0005524")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("purple", "GO:0005524")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("yellow", "GO:0046872")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("red", "GO:0019252")) -> MA_p_vals
c(MA_p_vals, pValReportingMA("yellow", "GO:0020037")) -> MA_p_vals

MA_enrichment_sizes %>% filter(n >= 50) %>% mutate("p-vals" = MA_p_vals) -> MA_enrichment_with_pvals
#read_tsv("E:/AshtonData/MPI_all_data/Datasets/GO_annots_descriptions.tsv") -> Go_anots

#inner_join(MA_enrichment_with_pvals, Go_anots, by = c("GO_ID")) -> ma_complete
#inner_join(rna_enrichment_with_pvals, Go_anots, by = c("GO_ID")) -> rna_complete

write_tsv(MA_enrichment_with_pvals,"E:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_GO_enrichment_p-vals.tsv")
write_tsv(rna_enrichment_with_pvals,"E:/AshtonData/MPI_all_data/GO_enrichment_analysis/RNA_GO_enrichment_p-vals.tsv")

