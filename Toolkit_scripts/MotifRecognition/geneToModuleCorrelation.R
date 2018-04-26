#!/usr/bin/env Rscript
#A tool to compare the correlation of a gene and the correlation of genes in a given gene module.
#switch key:

library("optparse")

##Parse input arguments
option_list = list(
  make_option(c("-g", "--genes"), default = NULL, help = "File reference to gene list you wish to check correlation of"),
  make_option(c("-m", "--module", default= NULL, help = "Path to the file containing the module genes")))
##Add an option to add a list of genes, instead of just one to speed up runing
opt_parser = OptionParser(option_list= option_list)
opt = parse_args(opt_parser)
if (is.null(opt$genes) || is.null(opt$module)){
  print_help(opt_parser)
  stop("Please supply valid arguments", call.=FALSE)
}


load(file = "/home/likewise-open/ICE/aomdahl/Datasets/RNASeq/CorrelationSuperMatrix.RData")
#write(paste("Loading genes from", opt$genes), file="")
lookup.genes.table = read.table(opt$genes, header = FALSE, stringsAsFactors=FALSE)
lookup.genes = lookup.genes.table$V1

module.genes = read.table(opt$module, header = FALSE, stringsAsFactors=FALSE)
#module.genes = read.table("/home/likewise-open/ICE/aomdahl/CoExpressionNetworks/MotifRecognition/TC-NIAT04811/NIATv7_g04811_module.txt", header = TRUE)

module.gene.correlations <- G_similarityF[module.genes$V1,]
#Grab the genes we want
average.module.gene.correlations <- as.vector(colMeans(module.gene.correlations))

correlationResults = vector()
#write.table(lookup.genes, file="")
for (gene in lookup.genes)
{
  #write(paste("Looking up gene", gene), file="")
  gene.vector <- as.vector(unlist(G_similarityF[gene, ]))
  #gene.vector = as.vector(unlist(G_similarityF[lookup.genes,]))
  final.cor <- cor(average.module.gene.correlations, gene.vector)
  #write(paste("Correlation calculated is", final.cor), file="")
  correlationResults[gene] <- final.cor
}

write.table(correlationResults, file= "", append = FALSE, row.names = TRUE, col.names=F, sep='\t', quote=FALSE)
write.table(correlationResults, file= paste(opt$module, "_correlation.tsv", sep=""), append = TRUE, row.names = TRUE, col.names=F, sep='\t', quote=FALSE)
