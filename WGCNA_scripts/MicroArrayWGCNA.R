#!/usr/bin/env Rscript
######################################################################
##  THIS VERSION OF THE SCRIPT OMITS COMMAND LINE OPTIONS
##  AND USES EXPRESSION DATA AND GINI CORRELATION
## 
##  Made by Ashton Omdahl, ICE, June-August 2016


#Simple tool that draws a density graph of points above
#A certian threshold

suppressMessages(library(WGCNA))
suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(DESeq2))
options(stringsAsFactors=FALSE);
enableWGCNAThreads()

STD_CORRELATION_THRESHOLD = 0.90
SOFT_POWER_THRESHOLD = 0.9
EXCEEDS_DISTANCE =0.15



##############Pre-processing of MicroArray Expression Data########################
ma.expression <- read.table("/home/likewise-open/ICE/aomdahl/Datasets/MicroArray/microArray_kinetics.tsv", header = TRUE, sep='\t')
#row.names(ma.expression) <- ma.expression[,"probeID"]
ma.datExpr <- as.data.frame(t(ma.expression))

#Check for genes/samples with too many missing values (should be good)
gsg = goodSamplesGenes(ma.datExpr, verbose=3)
gsg$allOK
##GOOD TO GO! Aug 9

#Do a simple hierarchical cluster and remove outliers from that test.
exprTree = hclust(dist(ma.datExpr), method = "average")
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(exprTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h=425, col="red")
#This looks pretty good, no visible outliers.- Aug 9

#programmatically, we want this to be anything outside of what, 2 sds??  need to adjust this.
#removingOutliers = cutreeStatic(exprTree, cutHeight = 425, minSize=5)
#keepSamples = (removingOutliers == 1)
#ma.datExpr <-datExprStabilized[keepSamples,]
#None removed!!

#Alternatively- this works, but always picks the furtherst.  not very useful.
#datExprFINAL2 <- datExprStabilized[!outlier(exprTree$height, logical=TRUE)]

save(ma.datExpr, file="Normalized_Expression_data.RData")


#MicroArray- which genes are most co-expressed
ma.datExpr.numeric <- data.matrix(ma.datExpr)


################################Identify most correlated genes ##################################################

#Read in correaltion data for reference
  #This is a 5gb file.  Keep only the stuff from here that we want- top 10,000
ma.G_similarity <- fread("/home/likewise-open/ICE/aomdahl/Datasets/MicroArray/kineticsCorrelationGCC.tsv", header = TRUE, sep='\t',col.names=T,data.table=F)
names(ma.G_similarity) <- row.names(ma.G_similarity)
ma.gcc_sim<-data.frame(G_similarity) 
row.names(ma.gcc_sim)<- names(ma.gcc_sim)
#Limit the genes shown in the expression data to those we are interested in and relevant to correlation data
ma.genes.shared<-intersect(names(ma.datExpr),row.names(ma.gcc_sim));

#These objects contain our key filtered data- highest correlations and stuff
ma.GCC.selection<-ma.gcc_sim[ma.genes.shared,ma.genes.shared]
ma.datExpr.selection <- ma.datExpr[,ma.genes.shared]

#Select the top 10,000 most interconnected genes from this overlapping subset
#Soft power threshold for us here- 0.9
maex.powers = c(c(1:10), seq(from =12, to=20, by=2))
ma.sft = pickSoftThreshold(ma.datExpr.selection, powerVector = maex.powers, RsquaredCut = SOFT_POWER_THRESHOLD)

ma.topGenes <-softConnectivity(ma.datExpr.selection, type="unsigned", power = ma.sft$powerEstimate, blockSize = 10000)

ma.topGenesWithNames <-as.data.frame(ma.topGenes)
rownames(ma.topGenesWithNames) <- names(ma.datExpr.selection) 

#Order them from most to least connected
ma.topTen = ma.topGenesWithNames[order(-ma.topGenes), , drop=FALSE]
#Keep only the top 10,000
ma.topTen = ma.topTen[1:10000,,drop=FALSE]
ma.topTen["NIATv7_g20693",] == ma.topGenesWithNames["NIATv7_g20693",]

#So the mapping is all good and solid- Aug 9
save(ma.topTen, file="TopTenThousandGenes.RData")
#########################################################################################


####################Isolate the top 10,0000 from the similarity matrix###################

ma.top.corr<- ma.GCC.selection[row.names(ma.topTen),row.names(ma.topTen)]

save(ma.top.corr, file = "TopTenThousandSimilarity_GCC.RData")
#save(ma.G_similarity, file="CorrelationSuperMatrix.RData")

#########Get the soft threshold based on this limited number of expressed genes ############
topCorrGenesNumeric <- matrix(as.numeric(unlist(ma.top.corr),rownames.force=TRUE),nrow=nrow(ma.top.corr))
ma.topCorrGenesNum <- topCorrGenesNumeric
g_powers = c(c(1:10), seq(from =12, to=30, by=2))
ma.g_sft = pickSoftThreshold.fromSimilarity(topCorrGenesNumeric, powerVector = g_powers, RsquaredCut = SOFT_POWER_THRESHOLD)
ma.g_softPower = as.numeric(g_sft$powerEstimate)
save(ma.g_sft, ma.g_softPower, topCorrGenesNumeric, file="MA_PreAdjacencyCalcs.RData")

#####Do the clustering################
#Create the adjacency matrix
ma.corr.adjacency <- adjacency.fromSimilarity(ma.topCorrGenesNum, type = "signed", power = ma.g_softPower)

#Topological overlap matrix
ma.g_TOM = TOMsimilarity(ma.corr.adjacency, TOMType="signed");
ma.dissTOM = 1-ma.g_TOM

#Cluster the differences in topological overlap to identify groups
ma.geneTree = hclust(as.dist(ma.dissTOM), method="average")
ma.geneTree$labels = (row.names(ma.topTen))

#Check where we are.
plot(ma.geneTree, xlab="", sub="", main = "MicroArray Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
######Checkpoint passed 8/10
###########Identify the modules#################

#Select the minimum size of modules- set to 50 to try and match the number in RNA-seq dataset
minModSize = 50

#Create the modules-- set deep split to 2 here to try and get closer to the number of modules we had in the other set.
ma.dynamicMods = cutreeDynamic(dendro=ma.geneTree, distM=ma.dissTOM,
                            deepSplit = 3, pamRespectsDendro=FALSE,
                            minClusterSize = minModSize)
#See the different sets:
table(ma.dynamicMods)

#Assign each gene co-expression module a color
ma.colors = labels2colors(ma.dynamicMods)

####Plot the modules###########
plotDendroAndColors(ma.geneTree, ma.colors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Modules:MicroArray")

save(ma.corr.adjacency, ma.g_TOM, ma.dissTOM, ma.geneTree, ma.dynamicMods, file="MicroArrayModules.RData")

#Passed checkpoint 8/10
############Record modules for output#############################
#Calculate intramodular Connectivity
ma.modularConn.scaled<-intramodularConnectivity(ma.corr.adjacency, ma.colors, scaleByMax= TRUE)

ma.modularConn.unscaled<-intramodularConnectivity(ma.corr.adjacency, ma.colors, scaleByMax= FALSE)

ma.modules.masterlist <- data.frame(Genes = c(row.names(ma.topTen)), group = c(ma.colors), 
                             ScaledConnectivity=ma.modularConn.scaled$kWithin, RawConnectivity=ma.modularConn.unscaled$kWithin)
row.names(ma.modules.masterlist) <- ma.modules.masterlist[,1]
ma.modules.masterlist<-ma.modules.masterlist[,-1]

ma.modules.ordered <- ma.modules.masterlist[rev(order(ma.modules.masterlist$group,ma.modules.masterlist$ScaledConnectivity)),]
fileName = paste("MicroArray_modules.tsv")
write.table(ma.modules.ordered, file= fileName, row.names = TRUE, sep = '\t')

##Checkpoint 8/10
#############################################################################################################
#Identify connectivity of TFs within Modules
#############################################################################################################

#Load the TF file
tf.genes <- read.csv("/home/likewise-open/ICE/aomdahl/Datasets/TranscriptionFactors/TranscriptionFactorFamilies.csv")
tf.genes <- as.data.frame(tf.genes)

#Keep only the TFs that are within modules
ma.tfs.in.modules<-intersect(row.names(ma.modules.ordered),tf.genes$Gene)

for (TF.gene.id in ma.tfs.in.modules){
  
  #Get the gene's group-- test case: NIATv7_g02976
  ma.module<-ma.modules.ordered[TF.gene.id,"group"]
  #select the gene's module- size 841, this is correct
  ma.modules.ordered.selected<-row.names(ma.modules.ordered[ma.modules.ordered$group==ma.module,]);
  #Get its rank from the ordered list
  sc.score <-as.numeric(ma.modules.ordered[TF.gene.id,][2])
  #Is the gene in the upper 50% of connectivity within the module?
  #if (rank <= length(modules.ordered.selected)/2)
  
  #Modification: is the scaled connectivity score > .5?
    if(sc.score > 0.5)
  {
    #Call gini lookup
    ## check connectivity;
    correlated.Genes <- ma.top.corr[TF.gene.id, ma.modules.ordered.selected]
    #Changed the size of the groups to 10%
    top.correlated.modular.genes <- (rev(sort(correlated.Genes)))[1:(0.1*length(correlated.Genes))]
    
    #If the top 10% doesn't give us a large enough list, take the top 30- I actually think 50 would be better (Aug 10)....
    if (length(top.correlated.modular.genes) <= 50)
    {
      top.correlated.modular.genes <- (rev(sort(correlated.Genes)))[1:50]
    }
    
      
    #Print this to an output file.
    fileName = paste("MicroArrayModules/",TF.gene.id,"_",ma.module, "_module.txt", sep="")

    focus.genes <- data.frame(Gene = c(names(top.correlated.modular.genes)), Correlation = as.numeric(top.correlated.modular.genes), 
                              ScaledConnectivity=c(ma.modules.ordered[names(top.correlated.modular.genes),]$ScaledConnectivity))
    #create the print table
    
    write.table(focus.genes, file= fileName, append = TRUE, row.names = FALSE, sep='\t', quote=FALSE, col.names=F)
    
    
  }
  else
  {
    next
  }
  
}

#############################Looking at specific tissue types

#Extract just the leaf data
#leaf.datExpr <-ma.datExpr.numeric[grep("SL|TL", row.names(ma.datExpr.numeric)),]
leaf.datExpr <-ma.datExpr.numeric[grep("SL|TL", row.names(ma.datExpr.numeric)),]  #Ran through on this one, August 11
root.datExpr <-ma.datExpr.numeric[grep("R", row.names(ma.datExpr.numeric)),]
#Do a simple hierarchical cluster and remove outliers from that test.
#Leaf data loks good, root data has an outlier.  Removed 8/10
testTree = hclust(dist(root.datExpr), method = "average")
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(testTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h=425, col="red")
removingOutliers = cutreeStatic(testTree, cutHeight = 80, minSize=5)
keepSamples = (removingOutliers == 1)
root.datExpr <-root.datExpr[keepSamples,]







##Doing this the fast way, August 10
tissue.powers = c(c(1:10), seq(from = 12, to=20, by=2))
leaf.sft = sft = pickSoftThreshold(leaf.datExpr, powerVector = tissue.powers)

root.sft = sft = pickSoftThreshold(root.datExpr, powerVector = tissue.powers)

leaf.adjacency = adjacency(leaf.datExpr, power = leaf.sft$powerEstimate)

root.adjacency = adjacency(root.datExpr, power = root.sft$powerEstimate)

leaf.TOM = TOMsimilarity(leaf.adjacency)
root.TOM = TOMsimilarity(root.adjacency)

leaf.dissTOM = 1-leaf.TOM
root.dissTOM = 1-root.TOM

leaf.geneTree = hclust(as.dist(leaf.dissTOM), method = "average")
root.geneTree = hclust(as.dist(root.dissTOM), method = "average")

tissue.minModuleSize = 50
leaf.mods = cutreeDynamic(dendro = leaf.geneTree, distM = leaf.dissTOM, deepSplit = 1, pamRespectsDendro= FALSE, minClusterSize = tissue.minModuleSize)
root.mods = cutreeDynamic(dendro = root.geneTree, distM = root.dissTOM, deepSplit = 1, pamRespectsDendro= FALSE, minClusterSize = tissue.minModuleSize)
#********
table(leaf.mods)
table(root.mods)

leaf.colors = labels2colors(leaf.mods)
root.colors = labels2colors(root.mods)
plotDendroAndColors(leaf.geneTree, leaf.colors, "Dynamic Tree Cut", dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05,
                    main = "Leaf gene module dendrogram")
plotDendroAndColors(root.geneTree, root.colors, "Dynamic Tree Cut", dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05,
                    main = "Root gene module dendrogram")

leaf.MEList = moduleEigengenes(leaf.datExpr, colors = leaf.colors)
root.MEList = moduleEigengenes(root.datExpr, colors = root.colors)
leaf.MEDiss = 1-cor(leaf.MEList$eigengenes)
root.MEDiss = 1-cor(root.MEList$eigengenes)
#****
leaf.METree = hclust(as.dist(leaf.MEDiss), method="average")
root.METree = hclust(as.dist(root.MEDiss), method="average")

plot(leaf.METree, main = "Clusting of eigengenes, leaf", xlab = "", sub = "")
#Very small height cut- maybe 0.1 or 0.05
sizeGrWindow(7,6)
plot(root.METree, main = "Clusting of eigengenes, root", xlab = "", sub = "")
#Also Very small height cut- maybe 0.1
tissue.METhresh = 0.15
abline(h=tissue.METhresh, col = 'red')


##Merge them
leaf.merge = mergeCloseModules(leaf.datExpr, leaf.colors, cutHeight = tissue.METhresh)
root.merge = mergeCloseModules(root.datExpr, root.colors, cutHeight = tissue.METhresh)
#Check out our work:
plotDendroAndColors(leaf.geneTree, cbind(leaf.colors,leaf.merge$colors), "Merged Tree cut", dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05,
                    main = "Leaf gene module dendrogram")
plotDendroAndColors(root.geneTree, cbind(root.colors,root.merge$colors), "Merged Tree cut", dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05,
                    main = "Leaf gene module dendrogram")
save(leaf.merge, root.merge, leaf.TOM, root.TOM,leaf.datExpr, root.datExpr, file = "TissueCoreData.RData")
################Calculating Intramodular connectivity################
leaf.conn.scaled<-intramodularConnectivity(leaf.adjacency, leaf.merge$colors, scaleByMax= TRUE)
root.conn.scaled<-intramodularConnectivity(root.adjacency, root.merge$colors, scaleByMax= TRUE)
leaf.conn.unscaled<-intramodularConnectivity(leaf.adjacency, leaf.merge$colors, scaleByMax= F)
root.conn.unscaled<-intramodularConnectivity(root.adjacency, root.merge$colors, scaleByMax= F)
#**

leaf.modules.masterList <- data.frame(Genes = c(colnames(leaf.datExpr)), group = c(leaf.merge$colors), 
                                    ScaledConnectivity=leaf.conn.scaled$kWithin, RawConnectivity=leaf.conn.unscaled$kWithin)
root.modules.masterList <- data.frame(Genes = c(colnames(root.datExpr)), group = c(root.merge$colors), 
                                      ScaledConnectivity=root.conn.scaled$kWithin, RawConnectivity=root.conn.unscaled$kWithin)


leaf.modules.masterList <- leaf.modules.masterList[rev(order(leaf.modules.masterList$group, leaf.modules.masterList$ScaledConnectivity)),]
root.modules.masterList <- root.modules.masterList[rev(order(root.modules.masterList$group, root.modules.masterList$ScaledConnectivity)),]

row.names(leaf.modules.masterList) <- NULL
row.names(root.modules.masterList) <- NULL

write.table(leaf.modules.masterList, file= "LeafModuleList.tsv", row.names = F, sep = '\t')
write.table(root.modules.masterList, file= "RootModuleList.tsv", row.names = F, sep = '\t')

leaf.corMatrix <- cor(leaf.datExpr)
root.corMatrix <-cor(root.datExpr)

############Now, actually creating the modules for use!
#Load the TF file
tf.genes <- read.csv("/home/likewise-open/ICE/aomdahl/Datasets/TranscriptionFactors/TranscriptionFactorFamilies.csv")
tf.genes <- as.data.frame(tf.genes)

#Alright, testing this

moduleBuilder <-function(orderedMList, corMatrixIn, tissueType){
  #Keep only the TFs that are within modules
  tfs.in.modules<-intersect(orderedMList$Genes,tf.genes$Gene)
  for (TF.gene.id in tfs.in.modules){
    print(TF.gene.id)
    #Get the gene's group-- test case: NIATv7_g02976
    curr.module<-orderedMList[orderedMList$Genes ==TF.gene.id,]
    #select the gene's module- size 841, this is correct
    selected.genes<-(orderedMList[orderedMList$group==curr.module$group,])$Genes
    corMatrix<- corMatrixIn[selected.genes,selected.genes]
    
  
    correlated.Genes <- corMatrix[TF.gene.id, selected.genes]
  
    #Changed the size of the groups to 10%
    top.correlated.modular.genes <- (rev(sort(correlated.Genes)))[1:(0.1*length(correlated.Genes))]
      
      #If the top 10% doesn't give us a large enough list, take the top 30- I actually think 50 would be better (Aug 10)....
      if (length(top.correlated.modular.genes) <= 50)
      {
        top.correlated.modular.genes <- (rev(sort(correlated.Genes)))[1:50]
      }
      
      
      #Print this to an output file.
      fileName = paste(tissueType,"/",TF.gene.id,"_",ma.module, "_module.txt", sep="")
      row.names(orderedMList) <- orderedMList$Genes
      focus.genes <- data.frame(Gene = c(names(top.correlated.modular.genes)), Correlation = as.numeric(top.correlated.modular.genes), 
                                ScaledConnectivity=c(orderedMList[names(top.correlated.modular.genes),]$ScaledConnectivity))
      #create the print table
      
      write.table(focus.genes, file= fileName, append = TRUE, row.names = FALSE, sep='\t', quote=FALSE, col.names=F)
      
      
    }
  
    
  }



moduleBuilder(root.modules.masterList, root.corMatrix, "RootTissueModules")
moduleBuilder(leaf.modules.masterList, leaf.corMatrix, "LeafTissueModules")