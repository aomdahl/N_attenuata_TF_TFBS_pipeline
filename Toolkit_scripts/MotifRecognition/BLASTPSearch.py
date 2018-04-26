#!/usr/bin/env python
#Formatting for SLimSuite CompariMotif tool
#
#Copyright 2016 aomdahl <aomdahl@ECO-12>
#

import sys
import os
import re
import argparse
import subprocess
import cPickle as pickle
import MotifList #??
import refMotif ##? from refMotif import refMotif.
import errorReporting as ER
sys.path.insert(0,"/home/likewise-open/ICE/aomdahl/Toolkit")
import HomerTools as HT
import GeneUtils as GU
	
####################################################################################################################


#Map of species names to abbreviations in the database for easy access
#Stored as {SpeciesName : [GeneHeader, gene abbreviation]
#TODO: Validate all these
speciesMap = {"Cucumis sativus" : ["Cucsa", "csa"], "Arabidopsis thaliana" : ["AT", "ara"], "Capsicum annuum": ["Capana", "can"],
              "Solanum melongena": ["Sme","sme"], "Solanum lycopersicum": ["Solyc", "sly"], "Nicotiana obtusifolia":["NIOBT", "nio"],
              "Vitis vinifera" : ["GSVIVT", "vvi"], "Glycine max" : ["Glyma", "gma"], "Mimulus guttatus" : ["Migut", "mgu"], "Nicotiana attenuata" : ["NIAT", "nat"],
              "Populus trichocarpa" : ["Potri", "ptr"], "Oryza sativa":["Os", "osa"]}

GENE_KEY = 0
ABBR = 1
global searchDB
searchDB = GU.BLASTP_PATH_NIAT

#List to regulate the reporting of a unavailable species to minimize output
globalReportingList = list()
#Simple method to check and see if a submitted
#GeneID and species are available in the database and match up
def checkGeneID(geneID, species, contRun = False):
	if not inSpeciesMap(species):
		if contRun:
			if species not in globalReportingList :
				print "BLASTp data unavailable for given species:", species
				globalReportingList.append(species)
			return False
		else:
			print "This species isn't available in the database."
			sys.exit()
	speciesName = inSpeciesMap(species)
	speciesGeneKey = speciesMap[speciesName][GENE_KEY]
	if speciesGeneKey not in geneID and speciesGeneKey not in unusalConversions(geneID):
		if contRun:
			#Write out the report to a file
			outStream = open("irregGeneIDBLAST.txt", 'a')
			outStream.write("Irregular GeneID encountered" +  geneID + ' '+ speciesName)
			#print "Irregular GeneID encountered", geneID, species
			return False
		print "This species doesn't seem to match the geneID you've given"
		print "Try the prefix id:", speciesGeneKey, "for", speciesName
		sys.exit()
	return True
	
	
	
#Checks to see if a species is actually in the species map or not.
#This requires its own function because its a bit more complicated- lists appear, etc.
#Returns the speices name or None if the name isn't there.
def inSpeciesMap(species):
	if species in speciesMap:
		return species
	sList = species.split(',')
	for s in sList:
		if s in speciesMap:
			return s.replace(" ", "")
	for stanSpec in speciesMap:
		for s in sList:
			if stanSpec in s:
				return stanSpec
	return None

#Makes the actual grep call
#@param geneID to search for
#@param species to specify the species, should make it go faster.
#@return a list of the hits
def BLASTPgrepSearch(geneID, species = None, db = GU.BLASTP_PATH_NIAT):
	try:
		#Clean up the geneID search query- sometimes, we get funny things....
		query = cleanQuery(geneID)
		if species:
			motifOverlaps = subprocess.check_output(["grep", species, db, "|", "grep", query]).split('\n')
		else:
			motifOverlaps = subprocess.check_output(["grep", query, db]).split('\n')

	except subprocess.CalledProcessError:
		return []
	return filter(None, motifOverlaps)


#A simple regular expression for checking geneIDs and cleaning them up, where possible
#@param geneID
#@return potentially cleaned up gene id.
def cleanQuery(geneID):
	cleanedID = unusalConversions(geneID)
	regex = "([\w\d\.\_]+)"
	reResponse = re.search(regex, cleanedID)
	if reResponse:
		return reResponse.group(1)
	else:
		print "Possible invalid geneID:", geneID
		return geneID


#This tool has a map of different gene nomenclatures
#and converts them to the format available in our BLASTp database
#2param the geneid
#@return the converted thing.
def unusalConversions(geneID):

	unusualConversions = {"testCase":"replace"}
	for word in unusualConversions:
		if word in geneID:
			print "Unusual case", word
			return geneID.replace(word, unusualConversions[word])
	return geneID


#This parses through and filters out the results of the BLASTp search
#@param results- grep results to actualy search
#@param responseNumber- the number of hits to return, script default is 1
#@param eThresh- the eThreshold to accept, default is all of them
#@param geneID- the gene ID we searched for
#@return the top responses in a list
def parseBLASTPResults(results, responseNumber, eThresh, geneID):
	if len(results) == 0:  #If we got no results, return empty table
		return results
		
	splitData = [x.split('\t') for x in results] #split up the grep hits by the columns
	splitData = [i for i in splitData if (geneID in i[0])]  #Keep those with the correct query search
	splitData.sort(key=lambda y:y[2], reverse = True)  #sort by the %identity
	if eThresh != -1:  #Keep only those passing the threshold as specified
		splitData = [z for z in splitData if float(z[10]) < float(eThresh)]
	#Returning the number back based on response number
	upperBounds = min(len(splitData)-1, responseNumber)
	return splitData[0:upperBounds]
	
	
#This function prints to the console the blast results with the following output:
#specifiedGene blastpHitGene %Similiarity %Value
#@param blastHits- the parsed grep search data.
#@return: All of the data printed to the console, if another function wants access to it.
#@return: this includes the [species, query, NIAT_response, %Overlap, e-value] 
#@return: If there are no hits, it returns none.
#@param printOut: if you actually want the data to print out, or just want it to return.
	
def printBLASTResults(blastHits, species, geneID,  printOut):
	retList =list()
	if not blastHits:
		if printOut:
			print "No results from search....", geneID
			open("NoBLASTPResults.txt", "a").write("No results from search: " +  geneID + '\n')
		#return None
		return [[species, geneID, "No hits", "", ""]]
	for hitData in blastHits:
		if printOut:
			print (species + '\t' + hitData[0][4:] + '\t' + hitData[1][4:]+'\t'+hitData[2]+'\t' + hitData[10])
		retList.append([species, hitData[0][4:], hitData[1][4:], hitData[2], hitData[10]])
	return retList
	
	
#This function makes the BLASTP database call, using a greative grep gethod :)
#@param geneID to search for 
#@param species we are interested in looking for
#@return a list containing all of the fianl responses.
def searchBLASTPDatabase(geneID, species, responseNumber, eThresh, contRun=False, printOut = True):
	#Make sure the species name is valid:
	species_f = species.replace("_", " ")
	query = cleanQuery(geneID)
	#Make sure the geneID is valid
	if not checkGeneID(query, species_f, contRun = contRun):  #Returns false if its a no good id.
		return
	#Run a grep search for that gene
	finalHits = parseBLASTPResults(BLASTPgrepSearch(query), responseNumber, eThresh, geneID)
	#Format this to a format that's pretty and nice for reading.
	return printBLASTResults(finalHits, species_f, geneID, printOut)
	
#This is another tool to optimize searching for lists of genes from the same species.
#Hopefully this bumps up our speed.
def searchListBLASTpBySpecies(species, geneList, responseNumber, eThresh, printOut = True):
	#Make sure the species names are valid
	species_f = species.replace("_", " ")
	retList = list()
	#Cont run must be true here to get through the entire list.
	cleanGeneList = [unusalConversions(i) for i in geneList if checkGeneID(i, species_f, contRun = True)]  #Keep those with the correct query search
	try:
		subSearchKey = speciesMap[inSpeciesMap(species_f)][ABBR]
	except KeyError:
		print "Species", species_f, "not in database"
		#Account for all those that aren't available.
		for i in cleanGeneList:
			reList.append([species_f, i, "Not in database", "", ""])
		return retList
		
#Preliminary grep call
	with open("blast_subsearch.txt", 'w') as subSearchFile:
		subprocess.call(["grep", subSearchKey, searchDB], stdout = subSearchFile)  #Create a subfile with the faster search material
		

	for geneID in cleanGeneList:
		query = cleanQuery(geneID)
 		finalHits = parseBLASTPResults(BLASTPgrepSearch(query, db="blast_subsearch.txt"), responseNumber,eThresh, query)
		tabularResults = printBLASTResults(finalHits, inSpeciesMap(species_f), geneID,  printOut)
		if tabularResults:
			retList.extend(tabularResults)
	#Remove the subsearch file- its quite big.
	subprocess.call(["rm", "blast_subsearch.txt"])
	return retList
	


if __name__ == '__main__':
	#Argument parsing
	parser = argparse.ArgumentParser(description='Performs a search of the NIAT BLASTp results.  If no option is specified, the program will automatically go from the local pickle dump file' + 
										          'Output prints to the console in format of: specifiedGene\tblastpHitGene\t%Similiarity\t%Value')
	parser.add_argument("-f", "--geneFile", help="Select this option to submit a file path with genes to search for their NIAT BLASTp match")
	parser.add_argument("-g", "--geneID",help = "Select this option to search for just a specific gene")
	parser.add_argument("-s", "--species", help = "The species from which the gene comes-- please give the full name: \"Genus species\"")
	parser.add_argument("-n", "--numberResponses",default = 1, type=int, help = "The maximum number of of blast hits to show.  Default shows only top hit.")
	parser.add_argument("-e", "--eThreshold", default = -1, type= float, help = "E-value threshold of values to accept, the default is all values.")
	parser.add_argument("--full", action = "store_true", help = "Specify if you wish to search the ALL x ALL BLASTp results.  Default is the NIAT hits only")
	parser.add_argument("-cp", "--consPVal", type = float, default = 1, help = "Specify a conservation p-value threshold for the motifs to search. Only for master list case. Default searches all")
	args = parser.parse_args()
	
	##Set the search database
	
	if args.full:
		searchDB = GU.BLASTP_PATH_FULL
		print "Using the full source database"
	
	
	if args.geneFile and args.species:
		globalReportingList = list()
		fileStream = open(args.geneFile, 'r')
		#parse the file
		#Call the response on each one
		#Spit out output
		print "This case"
	elif args.geneID and args.species:  #Only searching for a specific gene.
		#Search just for that one
		searchBLASTPDatabase(args.geneID, args.species, args.numberResponses, args.eThreshold)
		
	elif (args.geneFile or args.geneID) and not args.species:
		print "Please be sure to indicate the species from which the gene(s) come(s)."
		sys.exit()
	else:
		print "Default: using the serialzied data."
		masterMotifList = pickle.load(open(os.getcwd() + "/v3" + GU.SER_MOTIF_MAP, "rb"))
		masterMotifList.BLASTpSearchAllMotifs(args.numberResponses, args.eThreshold, consPValue= args.consPVal)
		GU.updateSerializedMotifData(os.getcwd() + "/vF" + GU.SER_MOTIF_MAP, masterMotifList)
	
	print "Unusual geneIDs written out to the file \'irregGeneIDBLAST.txt\'."
	
