library(readr)
library(dplyr)
library(tidyr)
library(magrittr)

rna_gene_list <- read_tsv("E:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_modules.tsv") %>% select("Gene")

ma_gene_list <- read_tsv("E:/AshtonData/MPI_all_data/GO_enrichment_analysis/ma_modules.tsv") %>% select("Gene")

go_reference_list <-read_tsv("E:/AshtonData/MPI_all_data/GO_enrichment_analysis/go_annots_cleanedup.tsv")

enrichment_permutation <- function(n, l_in)
{
  p_genes <- sample(l_in, n, replace = F)
  all_anots <- inner_join(data.frame("Gene" = p_genes), go_reference_list)
  all_anots %>% count(GO_ID) %>% arrange(desc(n))
}

enrichment_probabilities <- function(n, l_in, reps = 1000)
{
  final_t <- enrichment_permutation(n, l_in)
  for(i in 1:(reps-1)){
    final_t <-suppressWarnings(merge(x = final_t, y = enrichment_permutation(n, l_in), 
                                     by = c("GO_ID"), all = TRUE))
  }
  final_t[is.na(final_t)] <- 0
  final_t
}
sink("NUL")

#Added these in on 4/25, wanting to check the WRKY module in MA and the MYC2 modules out too
#57 in the MYC2 RNAseq
write_tsv(enrichment_probabilities(57, rna_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/MYC2_RNA_submodule_p_counts.txt")


#50 IN THE myc2 maSEQ
write_tsv(enrichment_probabilities(50, ma_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/MYC2_MA_submodule/MYC2_p_counts.txt")

#50 in teh WRKY MA
write_tsv(enrichment_probabilities(50, ma_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/WRKY3_MA/WRKY3_MA_submodule_p_counts.txt")





write_tsv(enrichment_probabilities(999, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_turquoise.txt")
write_tsv(enrichment_probabilities(982, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_blue.txt")
write_tsv(enrichment_probabilities(908, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_brown.txt")
write_tsv(enrichment_probabilities(846, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_yellow.txt")
write_tsv(enrichment_probabilities(784, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_green.txt")
write_tsv(enrichment_probabilities(749, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_red.txt")
write_tsv(enrichment_probabilities(667, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_black.txt")
write_tsv(enrichment_probabilities(604, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_pink.txt")
write_tsv(enrichment_probabilities(578, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_magenta.txt")
write_tsv(enrichment_probabilities(526, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_purple.txt")
write_tsv(enrichment_probabilities(473, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_greenyellow.txt")
write_tsv(enrichment_probabilities(279, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_tan.txt")
write_tsv(enrichment_probabilities(273, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_salmon.txt")
write_tsv(enrichment_probabilities(148, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_cyan.txt")
write_tsv(enrichment_probabilities(131, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_midnightblue.txt")
write_tsv(enrichment_probabilities(119, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_lightcyan.txt")
write_tsv(enrichment_probabilities(101, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_lightgreen.txt")
write_tsv(enrichment_probabilities(74, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_lightyellow.txt")
write_tsv(enrichment_probabilities(73, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_royalblue.txt")
write_tsv(enrichment_probabilities(72, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_darkred.txt")
write_tsv(enrichment_probabilities(63, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_darkgreen.txt")
write_tsv(enrichment_probabilities(58, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_darkturquoise.txt")
write_tsv(enrichment_probabilities(54, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_darkgrey.txt")
write_tsv(enrichment_probabilities(51, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_orange.txt")
write_tsv(enrichment_probabilities(47, rna_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/rna_darkorange.txt")
sink()


sink("NUL")
write_tsv(enrichment_probabilities(946, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_turquoise_.txt")
write_tsv(enrichment_probabilities(829, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_blue_.txt")
write_tsv(enrichment_probabilities(790, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_brown_.txt")
write_tsv(enrichment_probabilities(780, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_yellow_.txt")
write_tsv(enrichment_probabilities(677, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_green_.txt")
write_tsv(enrichment_probabilities(577, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_red_.txt")
write_tsv(enrichment_probabilities(552, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_black_.txt")
write_tsv(enrichment_probabilities(432, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_pink_.txt")
write_tsv(enrichment_probabilities(379, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_magenta_.txt")
write_tsv(enrichment_probabilities(370, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_purple_.txt")
write_tsv(enrichment_probabilities(362, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_greenyellow_.txt")
write_tsv(enrichment_probabilities(337, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_tan_.txt")
write_tsv(enrichment_probabilities(328, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_salmon_.txt")
write_tsv(enrichment_probabilities(245, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_cyan_.txt")
write_tsv(enrichment_probabilities(213, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_midnightblue_.txt")
write_tsv(enrichment_probabilities(206, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_lightcyan_.txt")
write_tsv(enrichment_probabilities(173, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_lightgreen_.txt")
write_tsv(enrichment_probabilities(156, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_lightyellow_.txt")
write_tsv(enrichment_probabilities(128, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_royalblue_.txt")
write_tsv(enrichment_probabilities(127, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_darkred_.txt")
write_tsv(enrichment_probabilities(119, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_darkgreen_.txt")
write_tsv(enrichment_probabilities(113, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_darkturquoise_.txt")
write_tsv(enrichment_probabilities(112, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_darkgrey_.txt")
write_tsv(enrichment_probabilities(110, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_orange_.txt")
write_tsv(enrichment_probabilities(106, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_darkorange_.txt")
write_tsv(enrichment_probabilities(106, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_white_.txt")
write_tsv(enrichment_probabilities(100, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_skyblue_.txt")
write_tsv(enrichment_probabilities(97, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_saddlebrown_.txt")
write_tsv(enrichment_probabilities(77, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_steelblue_.txt")
write_tsv(enrichment_probabilities(70, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_paleturquoise_.txt")
write_tsv(enrichment_probabilities(66, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_violet_.txt")
write_tsv(enrichment_probabilities(58, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_darkolivegreen_.txt")
write_tsv(enrichment_probabilities(53, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_darkmagenta_.txt")
write_tsv(enrichment_probabilities(206, ma_gene_list$Gene), "D:/AshtonData/MPI_all_data/GO_enrichment_analysis/MA_grey60_.txt")
sink()


#Now the code to run it for the submodules- on a whole scale. We should do on a miniscale too.
#Counts extracted by doing wc -l on module size.
write_tsv(enrichment_probabilities(30, rna_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g07696_darkgrey_module.txt")
write_tsv(enrichment_probabilities(52, rna_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g12711_purple_module.txt")
write_tsv(enrichment_probabilities(78, rna_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g16189_green_module.txt")
write_tsv(enrichment_probabilities(66, rna_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g17075_black_module.txt")
write_tsv(enrichment_probabilities(90, rna_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g18001_brown_module.txt")
write_tsv(enrichment_probabilities(57, rna_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g19088_magenta_module.txt")
write_tsv(enrichment_probabilities(57, rna_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g21131_magenta_module.txt")
write_tsv(enrichment_probabilities(57, rna_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g28241_magenta_module.txt")
write_tsv(enrichment_probabilities(78, rna_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g29978_green_module.txt")
write_tsv(enrichment_probabilities(47, rna_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g32572_greenyellow_module.txt")
write_tsv(enrichment_probabilities(90, rna_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g33798_brown_module.txt")
write_tsv(enrichment_probabilities(84, rna_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g34810_yellow_module.txt")

#MA stuff
write_tsv(enrichment_probabilities(78, ma_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g03410_yellow_module_ma.txt")
write_tsv(enrichment_probabilities(50, ma_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g04345_midnightblue_module_ma.txt")
write_tsv(enrichment_probabilities(50, ma_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g05680_magenta_module_ma.txt")
write_tsv(enrichment_probabilities(78, ma_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g15400_yellow_module_ma.txt")
write_tsv(enrichment_probabilities(50, ma_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g19262_magenta_module_ma.txt")
write_tsv(enrichment_probabilities(67, ma_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g19992_green_module_ma.txt")
write_tsv(enrichment_probabilities(55, ma_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g26487_black_module_ma.txt")
write_tsv(enrichment_probabilities(78, ma_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g27755_yellow_module_ma.txt")
write_tsv(enrichment_probabilities(82, ma_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g30725_blue_module_ma.txt")
write_tsv(enrichment_probabilities(50, ma_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g30738_royalblue_module_ma.txt")
write_tsv(enrichment_probabilities(78, ma_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g33088_yellow_module_ma.txt")
write_tsv(enrichment_probabilities(50, ma_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g40277_magenta_module_ma.txt")
write_tsv(enrichment_probabilities(50, ma_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g41088_magenta_module_ma.txt")
write_tsv(enrichment_probabilities(78, ma_gene_list$Gene), "E:/AshtonData/MPI_all_data/GO_enrichment_analysis/Submodule_analysis/NIATv7_g41243_yellow_module_ma.txt")

