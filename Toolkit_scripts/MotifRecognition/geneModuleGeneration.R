#!/usr/bin/env Rscript
#This tool generates a co-expression module for a given gene:

library("optparse")

##Parse input arguments
option_list <- list(
  make_option(c("-g", "--lookupGene"), help ="GeneID of the gene you wish to generate co-expression modules for. Omit the .t1 suffix"),
  make_option(c("-s", "--moduleSize"), default= 100, type="integer", help = "Size of the module you wish to extract, default is 100"),
  make_option(c("-n", "--name"), default = NULL, help = "Give the module a name if you want it"),
  make_option(c("-o","--output"), default= "./", help = "Output directory to write to.  Default is current working directory"),
  make_option(c ("-t", "--typeCluster"), default = "p", help="Specify which kind of module you want: pos. correlated genes (p), neg. correlated genes (n) or 1/2 of the top positives and 1/2 the top negatives (h)"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if(is.null(opt$lookupGene))
{
  print_help(opt_parser) 
}
if (is.null(opt$moduleSize)){
 
  print("Script will run with the following parameters:")
  print("Module size:")
  print(opt$moduleSize)
}
if(is.null(opt$name)){
  print("Module name")
  print(opt$lookupGene)
}

load(file = "/home/likewise-open/ICE/aomdahl/Datasets/RNASeq/CorrelationSuperMatrix.RData")

#Select the correlation values for the specified gene
module.gene.correlations <- G_similarityF[opt$lookupGene,]

#If there is no data there...
if (all(is.na(module.gene.correlations)))
{
  print("Unfortunately, data isn't available for that specified gene.")
  stop("Please make sure the gene name is correct", call.=FALSE)
}
gene.module <-vector()
switch(opt$typeCluster,
        
        n={gene.module = tail(sort(module.gene.correlations))[1:opt$moduleSize]},  #Negative options
        
        h={ upper <-tail(sort(module.gene.correlations, decreasing=TRUE))[1:(opt$moduleSize/2)]
           lower <-tail(sort(module.gene.correlations))[1:(opt$moduleSize/2)]
           gene.module = unlist(append(upper, lower))
        },
        
        p={gene.module = tail(sort(module.gene.correlations, decreasing=TRUE))[1:opt$moduleSize]},
        
        {gene.module = tail(sort(module.gene.correlations, decreasing=TRUE))[1:opt$moduleSize]}#Positives, the default
        ) 

#set the write out name.
moduleName <-  opt$name


if(is.null(opt$moduleName))
  {
  moduleName =opt$lookupGene
}

#Sort the model
gene.module <-sort(gene.module, decreasing=TRUE)
fileName = paste(opt$output,"/",moduleName, "_module.txt", sep="")
focus.genes <- data.frame(c(names(gene.module)), as.numeric(gene.module))
#create the print table
write.table(focus.genes, file= fileName, append = FALSE, row.names = FALSE, col.names=F, sep='\t', quote=FALSE)

print (paste("Gene module successfully written out to", fileName))
