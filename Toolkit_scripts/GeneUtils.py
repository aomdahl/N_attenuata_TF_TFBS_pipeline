#!/usr/bin/env python
###########################################
#This is a simple tool to perform basic, rendundant functions used throughout the motif prediction toolkit
##########################################
from Bio import SeqIO
import re
import gffutils
import glob
#import errorReporting as ER
import subprocess
import os
import cPickle as pickle
# A smple tool to flip a string



#Common keys for use throughout.
CS = "CandidateSequence"
DEATS = "Details"
EXACTS = "ExactMatchCount"
PVAL = "p-Value"
POS_FREQ = "PositiveFrequency"
NEG_FREQ= "NegativeBackgroundFrequency"
PMATRIX= "ProbabilityMatrix"
MOTIF_NOT_IN_OVERLAP = 0
MOTIF_IN_OVERLAP = 1
LODT = "LogOddsDetectionThreshold"
EVAL="e-value"
MATCH_SCORE = "MatchScore"
MATCH_RANK = "MatchRank"

ORTHO_INDEX = 0
OVERLAP_STATE_INDEX  = 1

NO_ORTHOLOG = "No Ortholog"
ORTHO_NO_OVERLAP = "Orthologous, but not conserved"
ORTHO_AND_OVERLAP = "Conserved"

EMPTY = ''
SER_MOTIF_MAP = "MOTIF_MAP.dump"
SER_GENE_MAP = "GENE_MAP.dump"
SUPER = "../"

PP ="PlantPan"
CISBP = "CISBP"
##############Database Paths
PP_TFs = "/home/likewise-open/ICE/aomdahl/Datasets/PlantPan/PlantPanMotifDB.tsv"
CISBP_TFs = "/home/likewise-open/ICE/aomdahl/Datasets/CIS-BP/CISBP_TF_DB.tsv"
GENOME_PATH= "/data/Genomes/NIATT/AssemblyNIATTr2/NIATTr2.fa"
ANNOTATION_PATH="/data/Genomes/NIATT/AssemblyNIATTr2/Annotation/Final/AN.6/NIATTr2.AN5_-_without_transposons.gff"
FUNC_ANNOTATION_PATH = "/data4/Results_recent/UPDATE_NEW_ANNOTATION/OUTPUT_BLAST/NEW_NIOBT_OUTPUT/blast2go/output/gene_annotation.NIATTr2_AN5.tsv"
BLASTP_PATH_NIAT = "/home/likewise-open/ICE/aomdahl/Datasets/ProteomeData/NIAT_subject_blastp.tsv"
BLASTP_PATH_FULL= "/data4/Results_recent/UPDATE_NEW_ANNOTATION/OUTPUT_BLAST/NEW_NIOBT_OUTPUT/FINAL_-_all-vs-all_-_ARA.CAN.CSA.GMA.MGU.NAT.NIO.PTR.SME.SLY.STU.VVI.tsv"
RNASEQ_MODULES = "/home/likewise-open/ICE/aomdahl/CoExpressionNetworks/GenomeWideAnalysis/TFBasedModules/GCC_modules.tsv"

ortho_DB_paths = {"Tomato":"/home/likewise-open/ICE/aomdahl/Datasets/GenomeData/NIAT_Orthologs/Tomato/NIAT-tomato-mapping.tsv", 
					"Arabidopsis":"/home/likewise-open/ICE/aomdahl/Datasets/GenomeData/NIAT_Orthologs/Arabidopsis/ARA-NATmapping.tsv", 
					"Poplar":"/home/likewise-open/ICE/aomdahl/Datasets/GenomeData/NIAT_Orthologs/Poplar/RBH-NAT.PTR-1e-6.tsv"} 

ortho_promoter_paths = {"Tomato":"/home/likewise-open/ICE/aomdahl/Datasets/GenomeData/NIAT_Orthologs/Tomato/SplitPromoterSequences/", 
					"Arabidopsis":"/home/likewise-open/ICE/aomdahl/Datasets/GenomeData/NIAT_Orthologs/Arabidopsis/SplitPromoterSequences/", 
					"Poplar":"/home/likewise-open/ICE/aomdahl/Datasets/GenomeData/NIAT_Orthologs/Poplar/SplitPromoterSequences/"} 

################Headers in the databse files####################
ID= "TF_ID"
MOTIF= "Motif_ID"
TF = "TF_Name"
SPECIES = "TF_Species"
FAMILY = "Family_Name"
DBDS = "DBDs"
DBID = "DBID"
###############Focus species list:
DEFAULT_SPECIES = ["Arabidopsis thaliana", "Cucumis sativus", "Populus trichocarpa", "Oryza sativa"]
FALSE_RESULTS = ["No", "No hits", "hits", "NA", "Data unavailable"]

def flipMinusStrand(strand):
    return strand[::-1]
    
def reverseSequence(seq):
	return flipMinusStrand(seq)

def genomeScaffoldMap(fileStream):
    genome = SeqIO.parse(fileStream, 'fasta')
    organismGeneMap={}
    for gene in genome:
        organismGeneMap[str(gene.id)] = str(gene.seq)
    return organismGeneMap



#Works in conjunction with Shuquing's python script
def extractPromoterSequence(fileIn):
    genome = SeqIO.parse(open(fileIn, 'r'), 'fasta')
    for gene in genome:
		if gene is None:
			continue
		print ">"+str(gene.id)
		promoterRegion = re.search("([a-z]*)[A-Z]", str(gene.seq))
		print promoterRegion.group(1).upper()
        
#Normalize score = (current score - minimum score) / (maximum score- minimum score)
#@param minScore
#@param maxScore
#@param score: the actual score
#@return the normalize score
def normalizeScore(minScore, maxScore, score):
	score_range = maxScore- minScore #The range between the score for exact match and minimum score, exact match score = float(query pattern length)
	relative_score = float(score) - minScore
	norm_score = relative_score/score_range
	return norm_score

#Uses a regular expression to extract a float percentage value from a string between parenthesis
def extractPercentage(stringIn):
	regexSearch = re.search('\((\d+\.\d+)%\)', stringIn)
	if regexSearch.group(1) is not None:
		return regexSearch.group(1)
	else:
		return "REGEX ERROR"


#This function identifies a p-value from a string based on its label and value
#At the moment, lazily made.  Make this more robust when you get a chance.
#stringIn the string with the P-value
#@return the p-value
def extractPValue(stringIn):
	pvalueS = stringIn[2:]
	return pvalueS


#Identifies the index of the column containing the desired trait
#matching score
#@param row: header
#@param word: to find in the header
#@return the index, or -1 if not match
def findSpecIndex(row, word):
	for x in range(0, len(row)):
		if row[x] == word:
			return x
			
	return -1

#A simple helper function to write out a list of items separated by tabs, with a newline at the end.
#@param fileStream to write to
#@param writeList to write out.
def writeTabList(fileStream, writeList):
	stringList = [str(x) for x in writeList]
	fileStream.write(('\t').join(stringList) + '\n')


def extractGeneAnnotations(nameList, outPath):
	out = open(outPath, 'w')
	db = gffutils.create_db(ANNOTATION_PATH,dbfn= "annotDB", force=True)
	for name in nameList:
		gene = db[name]
		print>>out, gene
		for i in db.children(gene, order_by='start'):
			print >> out, i
	out.close()

#This tool locates the local module file and returns a list of all the genes in the current module
#@param currentWorking diretory
#@return list5 of genes (strings)
def getModuleGenesList(searchDir, pathOnly = False):
	modulePaths = glob.glob(os.path.join(searchDir,'*_module.txt'))
	if len(modulePaths) == 0:
		##Look in the directory above
		print os.path.join(os.pardir,'*_module.txt')
		modulePaths = glob.glob(os.path.join(os.pardir,'*_module.txt'))
		if len(modulePaths) == 0:
			print "Unable to find module file- please ensure it is in the current or parent directory."
			return None
	if pathOnly:
		#print "Successfully found path:", str(modulePaths[0])
		return str(modulePaths[0])
	else:
		modulePath = modulePaths[0]
		fileStream = open(modulePath, 'r')
		geneList = list()
		for line in fileStream:
			table = line.strip().split('\t')
			if len(table) == 1:
				table = line.strip().split()
			if len(table) == 0:
				continue
			if table[0] != "Gene":
				geneList.append(table[0] + '.t1')
		return geneList


#This function gets the functional annotation of a gene from a gene annotation library
#@param geneName: to search for
#@return the annotation description
def getFunctionalAnnotations(geneName):
	geneID = removeGeneSuffix(geneName)
	try:
		search = subprocess.check_output(["grep", geneID, FUNC_ANNOTATION_PATH])
		search = search.strip().split('\t')
		if search is not  None: 
			if len(search) > 2:
				print "Multiple entries in geneName"
				return search[1] + search[2]
			return search[1]
	except subprocess.CalledProcessError:
		return ""


#Returns the module from which a gene comes
#@param geneName to search
#@param modulesDB- the module assignment db to examine, since there are multiple
#@param scaledConnectivity- if you want the scaled connectivity included as well.
#@return module assignment, none if none.
def getGeneModule(geneName, modulesDB = RNASEQ_MODULES, scaledConnectivity = False):
	geneID = removeGeneSuffix(geneName)
	try:
		search = subprocess.check_output(["grep", geneID, modulesDB]).strip().split('\t')
		if search is not None:
			if scaledConnectivity:
				try:
					return [search[1], search[2]]
				except IndexError:
					print "ScaledConnectivity not available for this module"
					return [search[1], "NA"]
			else:
				return search[1]
		else:
			return search[1]
			
	except subprocess.CalledProcessError:
		if scaledConnectivity:
			return ["No module assignment", ""]
		return "No module assignment"

#Method that removes the ".t1" if it is in the NIAt gene name
#@geneID to fix
#@return fixed gene
def removeGeneSuffix(geneID):
	if ".t1" in geneID:
		geneName = geneID.replace(".t1", "")
		return geneName
	else:
		return geneID

#def constructAnnotationsList(listFile, outPath):
def constructAnnotationsList():
	geneFile = open("./TestCase2/NIATTr2.FPKM.WWFAC.neg.txt", 'r')
	passList = list()
	
	for line in geneFile:
		line = line.strip()
		#print line
		if "NIAT" in line:
			phrase = re.search("NIAT\w*\s", line)
			passList.append((phrase.group(0))[:-1])
			print ((phrase.group(0))[:-1])
	extractGeneAnnotations(passList, "./TestCase2/negGenes.gff")
	geneFile.close()

#This is a general access function that calculates a conservation score (i.e. for a motif) as follows:

	#    Number of conserved pairs of the motif
	# -----------------------------------------------
	# Num. of ortholog pairs where the  motif appears
def calculateConservationScore(conservedPairs, orthologPairs):
	if orthologPairs != 0:
		return conservedPairs/float(orthologPairs)
	else:
		return -1
		print "No orthologous pairs found.."

#Tool to safely add to a map of keys to lists.
def safeAddtoListMap(mapObject, key_value, add_in):
	if key_value not in mapObject:
		mapObject[key_value] = list()
	mapObject[key_value].append(add_in)


##################################################GU toolkit for parsing databases.##################################################################################
#This tool searches the database for a certain term
#BRO: THIS TOOL IS GNARLY AWESOME

#This allows you to search a database for a term quickly
#@param dbPath: the database you want to search
#@param term: the term you are looking for(a string)
#@param searchKey: the header you want to search in for this term (a string)
#@param speciesFilter: a list of which species you want to include in the search, default is listed above
#@param returnKeys: which terms in the header you want back (a list)
#This can be easily optimized by only copying the data that actually has what you want instead of repeating.
#But this is a fast step anyway.
#@returns a map of the headers and their values.
def searchDatabase(dbPath, term, searchKey, returnKeys, speciesFilter = DEFAULT_SPECIES):
	#First, set up the species filters
	speciesFilterComplete = setUpSpeciesFilter(speciesFilter)
	db = open(dbPath, 'r')
	topLine = db.readline().strip() #get the top line to identify the headers
	db.close()
	headerMap = assignColumnKeys(topLine) # a map of columns to their indices

	#db = open(dbPath, 'r')
	##Modify this to use grep: more effecicnet:
	try:
		dbSearchResults = subprocess.check_output(["grep", term, dbPath]).split('\n')
	except subprocess.CalledProcessError:
		dbSearchResults = list()
	allResults =  lineSearchDatabase(term, headerMap[searchKey], filter(None, dbSearchResults), len(headerMap)) #Search all of the data for lines matching the term in the searchKey column
	#db.close()
	
	return filterResults(headerMap, allResults, speciesFilterComplete) #Filter the results just by these

#@param term to search for in...
#@param the searchKey column
#@param returnKeys the data to return
def searchPPDatabase(term, returnKeys, speciesFilter = DEFAULT_SPECIES):
	return searchDatabase(PP_TFs, term, MOTIF, returnKeys, speciesFilter)
	
def searchCISBPDatabase(term, returnKeys, speciesFilter = DEFAULT_SPECIES):
	return searchDatabase(CISBP_TFs, term, MOTIF, returnKeys, speciesFilter)
	
#This returns lines that match the term in the searchCol number
#@param term to search for
#@param searchCol to look at
#@param dbMatches- result of grep search on the database, faster than our linear search.
#@param length-number of columns
#@delim- delimter between data members
def lineSearchDatabase(term, searchCol, dbMatches, tableLength, delim = '\t'):
	if dbMatches is None or len(dbMatches) <= 0:
		return []
		
	retList = list()
	for line in dbMatches:
		data = line.strip().split(delim)
		#If the data is lacking stuff at the end for some reason,
		#need to append it so we don't get faults
		while len(data) < tableLength:
			data.append('')
		#print "Length fixed", str(len(data))			
		if data[searchCol] == term:
			retList.append(data)
	
	return retList


#returns a map of the column numbers to all the column headers
def assignColumnKeys(topLine, delim = '\t'):
	retMap = dict()
	tableTop = topLine.split(delim)
	for i in range(0, len(tableTop)):
		retMap[tableTop[i]] = i
	return retMap

#returnKeyIndices: the list of column numbers we want back.
#@return: a map of the headers and their values
#@param headerMap of the headers we want and their respective indices
#@param allResults is a list of all the matches
def filterResults(headerMap, allResults, speciesFilter):
	retList = list()
	for line in allResults:
		if speciesFilterCheck(line, speciesFilter):  #Does this incident contain our filtered species?
			subMap = dict()
			#Go through the mapping of headers : indices
			for key in headerMap:
				subMap[key] = (line[headerMap[key]])  #make the assignment
				#Special case for "details" column of CISBP --> replace with DNA binding domains (DBDS)
				if (key == DBDS):
					subMap[DEATS] = subMap[DBDS]
				#Special case for Gene id column in CISBP --> found in the DBID column.
				if (key == DBID):  
					subMap[ID] = subMap[DBID]

			retList.append(subMap)

	return retList
#A tool to check if a species is in our row to see if we want to keep it.
#@param line to check
#@param speciesFilter list of species looking for
#@return True if the match is there.
def speciesFilterCheck(line, speciesFilter):
		
	if any(y in line for y in speciesFilter):
		return True
		
	for phrase in line:
		for species in speciesFilter:
			if species in phrase:
				return True
				
	return False

#Converts the map of headers into a list of the columns we want.
def getColumnNumbers(returnKeys, headerKeys):
	retList = list()
	for element in returnKeys:
		if element in headerKeys:
			retList.append(headerKeys[element])
	return retList
#Helper function to parse out empty strings
def notEmptyString(strIn):
	if strIn != '' and strIn != "" and strIn is not None:
		return True


#Set up the species filter- DEFAULT_SPECIES needs additions, as does any custom input
#@parma the species list
#@returns the modified list.
def setUpSpeciesFilter(speciesList):
	listLen = len(speciesList)
	retList = list()
	retList.extend(speciesList)
	for term in speciesList:
		retList.append(term.replace(" ", "_"))
	return retList
	
#helper list to check a list, if its valid
#@param the list
#@return True if its good 
#@return false it its no good
def listChecker(listIn):
	if listIn is None:
		return False
	if len(listIn) == 0:
		return False
	if len(listIn) == 1 and listIn[0] in FALSE_RESULTS:
		return False
	return True

#A simple tool that looks up the correlation for a list of genes against a specified co-expression module
#@param geneList: the genes to lookup (typically putative TF genes)
#@param moduleFile: a file listing genes in the module
#@param db: specify the databse to use, standard is the RNA-seq one.  
def expressionCorrelationLookup(geneList, moduleFile, db = "RNA-seq"):
	#Put the genes into a file
	retMap = {}
	if not listChecker(geneList):
		return retMap
	tempFile = open("./putativeTFGeneList.txt", 'w')
	for gene in geneList:
		if gene not in FALSE_RESULTS:
			writeGene = removeGeneSuffix(gene)
			tempFile.write(writeGene + '\n')
	tempFile.close()
	
	#Call the R script
	try:
		mf = str(moduleFile)
		print "Examining module file", mf
		responseLines = subprocess.check_output(["Rscript", "/home/likewise-open/ICE/aomdahl/Toolkit/MotifRecognition/geneToModuleCorrelation.R", "-g","putativeTFGeneList.txt", "-m", mf])
		#print responseLines
		#raw_input("ahh.  this is why")
	except subprocess.CalledProcessError as err:
		print "Error in R script correlation extraction"
		print "Program will proceed with available results"
		responseLines = err.output
		
	responseCorrelations= responseLines.split('\n')
	for line in responseCorrelations:
		if "NIAT" in line:
			metrics = line.strip().split('\t')
			if len(line) > 1:
				retMap[metrics[0]] = metrics[1]
	return retMap
			
	#return a map of genes: correlation values(?)

#############################
#This is a useful tool for if we are running multiple runs of the test on a given module
#If updates the serialzed motif list by first checking if the list exists, and if it does, updating it
#Then, it writes out the file under the same name.
#Only works for the multi step approach (!)
#If you return things to one name, this has to go.
def updateSerializedMotifData(filePath, motifList, overwrite= False):
	#If the file has already been written, append it
	if not overwrite:
		if os.path.isfile(filePath):
			print "Merging current putative motif library with new results..."
			originalMotifList = pickle.load(open(filePath, "rb" ))
			motifList.append(originalMotifList)
		
	pickle.dump(motifList, open(filePath, 'wb'))
	print "Motif file library written out to", filePath
