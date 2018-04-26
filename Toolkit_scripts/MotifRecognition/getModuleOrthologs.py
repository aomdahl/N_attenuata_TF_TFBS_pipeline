#!/usr/bin/env python
#Formatting for SLimSuite CompariMotif tool
#
#Copyright 2016 aomdahl <aomdahl@ECO-12>
#
#

import sys
import os
import re
import argparse
import subprocess
import cPickle as pickle
from MotifList import MotifList
from enum import Enum
sys.path.insert(0,"/home/likewise-open/ICE/aomdahl/Toolkit")
import GeneUtils as GU
import random
import glob
import errorReporting as ER

#Class to help store the gene Orthologs and motifs
#Data structure basically associates a gene with its motifs
#And its ortholog
#
class geneOrtholog(object):
	
	def __init__(self, geneName):
		self.gene = geneName
		self.motifList = set()
		self.ortholog = None
		
	def addMotif(self, motif):
		#Makes sure the elements are unique!
		self.motifList.add(motif)

	def writeToFile(self, fileStream):
		fileStream.write(">" + self.gene + '\t'+'\t'.join(self.motifList))
		if self.ortholog is not None:
			fileStream.write("\nOrtholog: " + self.ortholog)
		fileStream.write('\n')
		
		
####################################################
#Key words
CNSRV_DIR = "./ConservationTesting"
INTERSECT = "intersectBed"
SIG_TESTING= "SignificanceTesting"
INTERSECT_CUTOFF ="0.9"

#Helper function to safely add a motif to the list
#@param dictionary of geneIDs to geneOrtholog objects, which contain its motif list and ortholog id.
#@param primary_key (here, the gene name)
#@param secondary_key for inner map (here, either MOTIF_SEQ or ORTHO)
#@param value_in is value to add

def addGeneMotif(map_in, gene_in, motif_in, updateSuper):
	if gene_in not in map_in:
		map_in[gene_in] = geneOrtholog(gene_in)
	map_in[gene_in].addMotif(motif_in)
	
	if updateSuper:
		#Keep track in the master list for later serialization.
		masterMotifList.getMotif(motif_in).addGene(gene_in)
	
#Simple tool using a rgular expression to extrac the gene name from the mapping.
#format of containing string is likely: three-key-abbreviation|geneID
#									ie: sly|Solyc05g047710.2.1
#@return the gene ID only
def getGeneName(containingString):
	regexSearch = re.search('\w\w\w\|([\w\.\d]+)', containingString)
	return regexSearch.group(1)


#This function searches the ortholog mapping object and finds the match.
#If a match is found, and the map is given, it adds the ortholog assignment to the gene
#@param geneMap -- map of all the genes associated with our module and its associated motifs
#@param geneKey- what gene we are searching for a mapping for.
#@param rand- do you just want a random ortholog?
#@return the ortholog name (string)
def findOrtholog(geneID, geneMap = None, rand = False):
	assert (orthologMapping is not None)
	#This will return a random Ortholog for our gene.  Happens in the case of our signficance testing.
	if rand is True:
		return random.choice(orthologMapping.values())
	else: 
		if (geneMap is not None) and (geneID in orthologMapping):  #Need to make sure its actually in there.
			setGeneOrtholog(geneID, geneMap,orthologMapping[geneID])
			return orthologMapping[geneID]
		return None


#Helper function to assign a gene object its ortholog
#Also updates the master list.
def setGeneOrtholog(geneID, geneMap, ortholog):
	#Add the ortholog to the key
	geneMap[geneID].ortholog = orthologMapping[geneID]
	#Update the master list
	masterMotifList.addOrtholog(geneID, orthologMapping[geneID])
	


#Loads the map into memory, storing it in a map object.
#Assumes the mapping is one-to-one.  If this is incorrect, then perhaps we have an issue.
def loadOrthologMap(orthoFilePath):
	orthologMap = {}
	orthoMappingFile = open(orthoFilePath, 'r')
	for line in orthoMappingFile:
		mappingLine = line.strip().split('\t')
		nGene = getGeneName(mappingLine[0])
		tGene = getGeneName(mappingLine[1])
		orthologMap[nGene] = tGene
	return orthologMap


#A tool to make a new directory safely.a
def makeNewDirectory(newDirect):
	#Set up new directory
	if os.path.isdir(newDirect):
		pass
	else:
		os.mkdir(newDirect)

#Gets out a motif from the .bed format..  Pretty simple.
def extractMotif(data_in):
	regexSearch = re.search('\d\-([A-Z]+)', data_in)
	if regexSearch.group(1) is None:
		return "NONE"
	return regexSearch.group(1)

#Simple function to loop through the items in the BED map 
#and write them to local files for quick access.
#@param geneKey: the current gene
#@param geneMotifBED: map of genes : BED data.
def writeBedData(geneKey, geneMotifBED):
	motifBEDFile = open(geneKey + "_motifs.bed", 'w')
	for motifLine in geneMotifBED[geneKey]:
		motifBEDFile.write(motifLine)
	motifBEDFile.close()

#This function parses the args.NIAT_genes input file, the output of scanMotifGenomeWide.pl -bed
#@return a list with 2 maps: [0] is a map of the gene and its motifs
#							 [1] is a map of each gene and the lines of text it has associated with it (used to split the .bed file)
#@param topMotifPath- path to the NIAT_genes input file.
def parseNIATMotifLocations(topMotifPath, updateSuper = True):
	
	motifGenesFile = open(topMotifPath, 'r') 
	motifGeneMap = {} #dictionary of geneIDs to geneOrtholog objects, which contain its motif list and ortholog id.
	
	lineMap = {} ##This stores the lines of the file for splitting at a later time.
	
	#Assign each gene its motifs
	for line in motifGenesFile:
		lineTable = line.strip().split('\t')
		gene = lineTable[0]
		
		#Keep track of the lines in the lineMap --MODIFY THIS TO USE VERSION IN GU
		if gene not in lineMap:
			lineMap[gene] = list()
		lineMap[gene].append(line)
		
		motif = extractMotif(lineTable[3])
		addGeneMotif(motifGeneMap, gene, motif, updateSuper)
	motifGenesFile.close()
	return [motifGeneMap, lineMap]

#CAlls the bedtools intersect function to identify intersection between 
#two files, and stores the output motifs in a map for further processing
#@param align_1: first file to align, here the gene/ortholog alignment file
#@param align_2, here the motifAlignment.bed file
#@param conservedMotifs- the map of conserved motifs thus far.
#@return map of Conserved motifs and their count frequency
def updateConservedMotifMap(align_1,  align_2, conservedMotifs, noOverlapsList = None, updateMasterList = True):
	MOTIF_NOT_IN_OVERLAP = 0
	MOTIF_IN_OVERLAP = 1
	#If motifOverlaps has a length of 0, there are no places where the motif overlaps the gene and its ortholog
	#These need to be stored in a list.
	motifOverlaps = intersect(align_1, align_2)
	#Assuming we got something back from the intersect call.
	if checkConservationOnAlignment(motifOverlaps, noOverlapsList):
		for line in motifOverlaps:
			intersectResults = line.strip().split('\t')
			if intersectResults[0] == GU.EMPTY:  #If we are at the last line, since there is an empty line.
				continue
			else:
				motif = extractMotif(intersectResults[3])
				gene = intersectResults[0]
				#Indicate that the motif is conserved
				if motif not in conservedMotifs:
					conservedMotifs[motif] = 0
				conservedMotifs[motif] += 1
				
				if updateMasterList:
					#Update this in the master list.
					masterMotifList.getMotif(motif).updateConservationExplicit(gene, MOTIF_IN_OVERLAP)
					
	return conservedMotifs

#CAlls the bedtools intersect function to identify intersection between 
#two files, and stores the output motifs in a map for further processing
#@param align_1: first file to align, here the gene/ortholog alignment file
#@param align_2, here the motifAlignment.bed file
#@return list of the output lines
def intersect(align_1,  align_2):
	#Arguments are given as follows:
	#Intersect: calls the bedtools intersect function
		#-b gene_ortho_alignment: the alignment bed file of the gene promoter and its ortholog promoter
		#-a motifs_alignment: the alignment bed file of the motifs in the gene promoters
		#-f minimum overlap required as a fraction of A (motifs_alignment file), must be nearly the whole motif.
	motifOverlaps = subprocess.check_output([INTERSECT,"-a",(align_2+ "_motifs.bed"), "-b",(align_1+ "_alignment.bed"), "-f", INTERSECT_CUTOFF]).split('\n')
	return motifOverlaps



#Boolean function that indicates if data from intersectBed have overlaps
#Updates the noMotifOverlaps list if there is no overlap at the location of the motifs.
#@param intersectOutput: output from Bedtools intersect function
#@param noMotifOverlaps: list of genes containing orthologs but no overlaps with motif and ortholog
#@return True if genes have intersection (i.e. the motif is conserved.
def checkConservationOnAlignment(intersectOutput, noMotifOverlaps):
	if intersectOutput[0] == GU.EMPTY:
		if noMotifOverlaps is not None:
			#print "Valid overlap list" ##
			noMotifOverlaps.append(currentGene)
		return False
	return True
	
#Makes the call to the command line to run the YASS program
#Prepares for this by creating a directory and copying the promoter files into it first
#Its pretty sweeeet.
#Genekey: the gene name we wish to align with...
#geneOrtholog: the ortholog we wish to align with
#orthologPromoterPath: path to access the ortholog's promoter
def runYASSAlignment(geneKey, geneOrtholog, NIATPromotersPath, orthologPromotersPath, folderName = None):
	gene_orthoName = geneKey + '_' + geneOrtholog  #A combined name for file writing
	if folderName is None:
		folderName = gene_orthoName
	makeNewDirectory(folderName)
	genePromoterNewFile = createLocalPromoterCopies(geneKey, NIATPromotersPath, folderName)
	orthologPromoterNewFile= createLocalPromoterCopies(geneOrtholog, orthologPromotersPath, folderName)
	os.chdir(folderName)
	if genePromoterNewFile is None or orthologPromoterNewFile is None:
		return None
	#print "yass ./" + genePromoterNewFile + " ./" + orthologPromoterNewFile + " -o " + gene_orthoName + "_alignment.bed -d 4"
	os.system("yass ./" + genePromoterNewFile + " ./" + orthologPromoterNewFile + " -o " + gene_orthoName + "_alignment.bed -d 4")
	return True

#Copies a file from a database of file split by gene [promoter sequence, .fsa files] to another directory [as a .fa file for YASS analysis]
#@param geneID: the geneID of the file to copy
#@param promoterDBPath: the path to the DB
#@param newDirectory: the new directory to write to
#@return the new file Name
def createLocalPromoterCopies(geneID, promoterDBPath, newDirectory):
	#Copy files to local directories
	promoterSourceFile = geneID + ".fsa"  #File name for the promoter original file
	promoterNewFile = geneID+ ".fa"   #File name for the local file copy
	#Check to see if its there already; if so, don't bother copying it.
	try:
		newDir = os.path.join(os.getcwd(), newDirectory)
		#raw_input(newDir)
		if not os.path.isfile(newDirectory + "/"+ promoterNewFile):
			subprocess.call(["cp", (promoterDBPath + '/' + promoterSourceFile), (newDir + "/"+ promoterNewFile)])
	except subprocess.CalledProcessError:
		return None
	except OSError:
		raw_input("error... :(")
		print "Unable to find the directories:"
		print "Promoter DB path", promoterDBPath
		print "Promoter Source File", promoterSourceFile
		print "New directory", newDirectory
		print "New promoter file copy", promoterNewFile
	return promoterNewFile

#Is there a better place for this?
def findSignificanceOfConservationScores(promoterDBPath, ortho_DB_path):
	makeNewDirectory(SIG_TESTING)
	os.chdir(SIG_TESTING)
	sigTestingDirectory = os.getcwd()
	motifConservationScores = {} #Storage place for the list of motif conservation scores
	
	for putativeMotif in masterMotifList.motifMap:  #Being done the easy way- try our custom iterator sometime!
		pMotifObject = masterMotifList.getMotif(putativeMotif)
		
		#orthologousGenes are genes that have orthologs
		orthologousGenes = pMotifObject.getGenesWithOrthologs()

		#Step through the list of orthologous genes and run the appropriate tests
		makeNewDirectory(sigTestingDirectory+"/" + putativeMotif)
		os.chdir(sigTestingDirectory+"/" + putativeMotif)

		comparisonDirectory = os.getcwd()
		
		#Run yass on each gene and its random ortholog 10 times
		for i in range(0,args.reps):  
			conservationMap = {putativeMotif:0}
			t_consCount = 0
			for gene in orthologousGenes:
				os.chdir(comparisonDirectory)
				randomOrtholog = findOrtholog(gene, rand = True)

				successful = runYASSAlignment(gene, randomOrtholog, promoterDBPath, ortho_DB_path, folderName = gene)
				if not successful:
					continue
				#Now, we want to run the intersction tool on this alignment and the motif's 
				#This requires a little extra work to get the path.
				gene_orthoName = gene + '_' + randomOrtholog  #A combined name for file writing
				
				originalMotifAlignmentPath = glob.glob(consPath + '/'+ gene +"_*")[0] + '/' + gene + "_motifs.bed"
				newMotifAlignmentPath = comparisonDirectory + "/" + gene + "/" + gene
				
				extractMotifAlignmentData(putativeMotif,originalMotifAlignmentPath, (newMotifAlignmentPath+ "_motifs.bed"))

				conservationMap = updateConservedMotifMap(gene_orthoName,newMotifAlignmentPath, conservationMap, updateMasterList = False)
	
			assert(len(conservationMap.keys()) <= 1)  #If this isn't true, we are looking at multiple motifs at a time.  uh oh.
			conservationScore = GU.calculateConservationScore(conservationMap[putativeMotif], len(orthologousGenes))
			GU.safeAddtoListMap(motifConservationScores, putativeMotif, conservationScore)  #For printing (?)
			masterMotifList.getMotif(putativeMotif).conservationScore
			
			
	#compare significance values here.
	#print motifConservationScores
	return motifConservationScores
	
	
#Selects lines from an alignment analysis pertaining only to a motif we are interested in.
def extractMotifAlignmentData(motif, in_path, out_path):
	#print "Original file:", in_path
	#print "New file:", out_path
	originalAlignment = open(in_path, 'r')
	filteredAlignment = open(out_path, 'w')
	for line in originalAlignment:
		line = line.strip()
		lineMotif = extractMotif(line.split('\t')[3])
		if lineMotif == motif:
			filteredAlignment.write(line + '\n')
	originalAlignment.close()
	filteredAlignment.close()

#Calculates the custom "p-value" for the conservation scores of each motif.  Updates the master motifList as well.
#@param bootstrapConservationScores: motifs and a list of bootstrap conservation scores
#@param masterMotifList: the master list with all the motif data.
#@return the p-value
def calculateConservationPVal(bootstrapConservationScores, masterMotifList):
	customP = -1
	for motif in bootstrapConservationScores:
		pMotifObject = masterMotifList.getMotif(motif)
		posConsScore = pMotifObject.conservationScore
		
		if posConsScore == -1:  #This means the motif isn't conserved; this call shouldn't be made at all.
			ER.storeError("Get Module Orthologs", "Conservation p-val calculation", "Unexpected result")
			return None
			
		#Gets the number of conservation scores greater than or equal to the positive data's conservation score.
		aboveConsScore = sum(i >= posConsScore for i in bootstrapConservationScores[motif])
		#What fraction of the bootstraps are above the positive conservation score?
		customP = aboveConsScore/float(len(bootstrapConservationScores[motif]))
		#If there were none:
		if customP == 0.0:
			customP = 1/float(len(bootstrapConservationScores[motif]))
			pMotifObject.setConsPValue(customP, "<")
		else:
			pMotifObject.setConsPValue(customP, "=")
	return customP


#This script generates directories in which YASS analysis can be run and then runs the YASS analysis
#Each directory has a sequence from the organism and the ortholog sequence from a comparison model organism.
global currentGene
global currentOrtholog
global noMotifOverlapsList
global orthologMapping
global consPath

if __name__ == '__main__':
	#Argument parsing
	parser = argparse.ArgumentParser(description='This script matches NIAT genes to their Tomato orthologs and extracts the corresponding promoter sequences from each.')
	parser.add_argument("NIAT_promoters", help="Path to the NIAT promoter sequence database, all written to individual files")
	parser.add_argument("NIAT_genes", help="Path to the .tsv file containing NIAT genes and their motifs.  Output of scanMotifGenomeWide.pl -bed")
	parser.add_argument("ortho_db", help="DB of ortholog's promoter sequences, all written to individual files.")
	parser.add_argument("i_map", help="One-to-one (injective) ortholog mapping for organism.")
	parser.add_argument("-s", "--sig_test", choices = ["y", "n"], default = "y", help = "Specify if you wish to perform conservation significance testing (y/n)")
	parser.add_argument("-r", "--reps", default = 1000, type=int, help = "number of times to repeat the randomized p-value testing.  Convenient here for quicker runs.")
	parser.add_argument("-p", "--pval", default = 1.0, type = float, help = "Specify the p-value threshold for writing out acceptable motifs.  Default is all")
	args = parser.parse_args()
	if not os.path.isdir(args.NIAT_promoters):
		print "Please ensure that the path to the NIAT promoters is in fact a directory, not a file"
		sys.exit()
	
	global masterMotifList
	#De-serialize motif data for access.
	masterMotifList = pickle.load(open("v1" + GU.SER_MOTIF_MAP, 'rb'))
	abspath = os.path.abspath("v1" + GU.SER_MOTIF_MAP)

	#Set up new directory
	makeNewDirectory(CNSRV_DIR)
	
	#Set up global variables and paths	
	noMotifOverlapsList = list()
	NIAT_genes_path = os.path.abspath(args.NIAT_genes)
	NIAT_promoters_path = os.path.abspath(args.NIAT_promoters)
	ortho_DB_path = os.path.abspath(args.ortho_db)
	
	#Load the ortholog map into memory
	mapPath = os.path.abspath(args.i_map)
	orthologMapping = loadOrthologMap(mapPath)
	
	#Create a map between each gene in the module and motifs it contains
	motifMaps = parseNIATMotifLocations(NIAT_genes_path)
	motifGeneMap = motifMaps[0] #this contains a list of all the motifs and their associated genes.
	geneMotifBED = motifMaps[1] #This is a map with a gene:its BED file contents for later use.
	
	#Track those genes which have no ortholog
	noOrtholog = list()
	#Keep track of 
	totalConservation= dict()

	#Move in to our directory where we will do the bulk of the work
	os.chdir(CNSRV_DIR)
	consPath = os.getcwd()
	
	#Loop through the list of genes
	for geneKey in motifGeneMap:
		currentGene = geneKey
		#First: find the orthologs for each gene, update the ortholog in to the map.
		findOrtholog(geneKey, geneMap = motifGeneMap)
		
		#If there is no ortholog, keep track of that in a separate list.
		if motifGeneMap[geneKey].ortholog is None:#This occurs if there is no ortholog.
			noOrtholog.append(motifGeneMap[geneKey])
			currentOrtholog = None
		else:  
			currentOrtholog = motifGeneMap[geneKey].ortholog
			successful = runYASSAlignment(geneKey, currentOrtholog, NIAT_promoters_path, ortho_DB_path)
			if not successful:  # if we are unable to run the yass
				continue
			#Write the BED data out here for intersection comparison
			writeBedData(geneKey, geneMotifBED)
			gene_orthoName = geneKey + '_' + motifGeneMap[geneKey].ortholog
			totalConservation = updateConservedMotifMap(gene_orthoName, geneKey, totalConservation, noOverlapsList = noMotifOverlapsList)
			os.chdir("../")
			
	masterMotifList.calculateConservationScores()
	#masterMotifList.writeToFile("./MOTIF_CONSERVATION_COMPLETE.txt")
	
	#Write genes with no orthologs to a file
	noOrthologFile = open("genes_no_ortholog.tsv", 'a')
	for gene in noOrtholog:
		gene.writeToFile(noOrthologFile)
	noOrthologFile.close()

	#Final Report
	print len(noOrtholog), " out of" , len(motifGeneMap), "genes in the module had no ortholog mapping."
	print len(noMotifOverlapsList), "genes had orthologs but had no intersection of a motif with the ortholog."
	print "For complete data results and conservation scores, see \"MOTIF_CONSERVATION_COMPLETE.txt\"."
	if(args.sig_test == "y"):
		bootstrapConservationScores = findSignificanceOfConservationScores(NIAT_promoters_path, ortho_DB_path)
		calculateConservationPVal(bootstrapConservationScores, masterMotifList)
	#Final file write outs	
	masterMotifList.writeToFile(consPath+"/MOTIF_CONSERVATION_COMPLETE.txt")
	#Write out a list of the motifs in homerFormat that make the cut.
	print "Writing out the conserved motifs to a .motifs file at ", consPath, "/conservedMotifs.motifs"
	masterMotifList.writeListHomerFormat(consPath + "/conservedMotifs.motifs", masterMotifList.getSigConserved(conservationP = args.pval))
	pickle.dump(motifGeneMap, open(GU.SER_GENE_MAP, 'wb'))
	#Keep track of all the data, so its available for later
	GU.updateSerializedMotifData(abspath.replace("v1", "v2"), masterMotifList, overwrite=True)
	#pickle.dump(masterMotifList, open(abspath.replace("v1", "v2"), 'wb'))

