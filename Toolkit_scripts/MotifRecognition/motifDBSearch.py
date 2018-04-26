#!/usr/bin/env python
#Formatting for SLimSuite CompariMotif tool
#
#Copyright 2016 aomdahl <aomdahl@ECO-12>
#

import sys
import os
import sys
import re
import argparse
import subprocess
import cPickle as pickle
from MotifList import MotifList
from refMotif import refMotif
import errorReporting as ER
sys.path.insert(0,"/home/likewise-open/ICE/aomdahl/Toolkit")
import HomerTools as HT
import GeneUtils as GU
	
#This function iterates through the directories in the comparison folder
#...looks at all the homer results, and creates refMotif objects for each one
#... and then stores a list of refMotif objects in a map associated with each motif
#@param path of HOMEr output files to search
#@param db to extract from
#@return a map of reference Motif objects.
def extractComparisonResults(path, db, speciesFilter):
	refMotifMap = dict()
	print "Writing out list of tfbs matches..."
	matrixMatches = open("known_tfbs_hits_" +db+ ".txt", 'w')
	for i in os.listdir(path):
		if i.endswith(".info.html"): 
			referenceMotifList = list()
			homologousTFBS = HT.motifMatchHTMLExtraction(path + "/" + i)
			motifSeq = homologousTFBS[0]
			referenceMotifs = homologousTFBS[1]  #The list contains [name, rank, score]
			referenceMotifs = filterTopResults(referenceMotifs, args.top)  #We actually want to keep them all until the final filtering, because most will be lost anyway.
			
			matrixMatches.write(motifSeq + '\n')
			for r in referenceMotifs:
				matrixMatches.write(r[0] + '\t' + r[1] + '\t' + r[2] + '\n')
				addInRef = DBLookUp(r, db,speciesFilter)
				if addInRef is not None:
					referenceMotifList.append(DBLookUp(r, db,speciesFilter)) 
			refMotifMap[motifSeq] = referenceMotifList
	matrixMatches.close()
	return refMotifMap
	


#Makes the database lookup call and stores the results in a refMotif object
#@param the reference motif name
#@param db to look in
#@return refMotif object
#TODO: is the scoping right here?  it doesn't look like it.
def DBLookUp(referenceMotif, db, filterSpecies):  #Rember, its [name, rank, score]
	refMotifList = list()
	if db == GU.PP:
		refMotifList = GU.searchPPDatabase(referenceMotif[0], PP_DEFAULT_RETURN_KEYS, speciesFilter = filterSpecies)
		
		if len(refMotifList) > 0:  
			refMotifObject = refMotif(refMotifList, GU.PP)
		
	else:
		refMotifList = GU.searchCISBPDatabase(referenceMotif[0],CISBP_DEFAULT_RETURN_KEYS, speciesFilter = filterSpecies)
		
		if len(refMotifList) > 0:
			
			refMotifObject = refMotif(refMotifList, GU.CISBP)

	if (len(refMotifList) <= 0):
		print "No match in the " + db +  " DB for orthologs of", referenceMotif[0] , "with the current plant filter selection."
		return None
		
	refMotifObject.setMatchParams(referenceMotif[1],referenceMotif[2])
	return refMotifObject


	
#This function keeps only the top motifs from the list, and makes sure they are sorted by score
#@param reference Motifs- the list to sort by match score
#@param topCount: the max number of motifs to keep 
#@return a list of the motifs with the highest match score, up to length topCount
def filterTopResults(referenceMotifs, topCount):
	sortedRefMotifs = sorted(referenceMotifs,key=lambda ref:ref[2], reverse = True) #Rember, its [name, rank, score]
	filteredRefMotifs = [i for i in referenceMotifs if float(i[2]) >=  float(args.matchScore)]

	if len(filteredRefMotifs) >= topCount:
		return filteredRefMotifs[:topCount]
	else:
		return filteredRefMotifs
		
#This only works if the maps have the same keys!
def combineReferenceMotifMaps(m1, m2):
	for motifSeq in m1:
		if motifSeq in m2:
			m1[motifSeq].extend(m2[motifSeq])

		else:
			print "Error, " , motifSeq, "..not in 2nd map"
		
		
	return m1

#Does what the name says: writes out the list of refernce motifs to a file.
#@param filePath : to write out to
#@param refMotifMap: a map of motif sequences : list of reference motifs
def writeRefMotifs(filePath, refMotifMap):
	print "Writing out to:", filePath
	outFile = open(filePath, 'w')
	#print refMotifMap
	for r in refMotifMap:
		outFile.write("####################################################\n")
		outFile.write(r + '\n')
		currList = refMotifMap[r]
		#print currList
		for refObj in currList:
			#print refObj
			if type(refObj) is list:
				print "We have a list object"
				print refObj
				print "Current list:"
				print currList
				raw_input("What's next?")
				
			refObj.writeToFile(outFile)
		
	outFile.close()
			

#Fairly simple: loops through the map of pMotifSequences : reference motifs
# and assigns each corresponding motif in the master list its new reference motifs
#@param masterList- the pickling master list
#@param refMotifList: the map of referenceMotifs
def updateMasterList(filterSpecies, masterList, refMotifMap):
	for pMotif in masterList.motifMap:
		if pMotif not in refMotifMap:
			print "It appears that the following motif is missing from the list of reference Motifs..."
			print pMotif
			print "This likely a result of HOMER merging similar motifs in the compareMotifs analysis."
			print "Analysis will continue.  To run further analysis on all motifs, increase the \' -reduceThresh\' variable"
			print "on the compareMotifs.pl call"

		else:
			masterList.motifMap[pMotif].refMotifs = refMotifMap[pMotif]
			#for debugging
#			for r in refMotifMap[pMotif]:
#				r.printOut()
	masterList.referenceSpecies = filterSpecies

	GU.updateSerializedMotifData("v3" + GU.SER_MOTIF_MAP, masterList)

#A tool to help set the ouput print file name, 
#@pparam dbName, the db you are in
def outFileName(dbName):
	PP_OUT = "./TF_matches_PP.txt"
	CISBP_OUT = "./TF_matches_CISBP.txt"
	DEFAULT = "./TF_matches.txt"
	if dbName == GU.PP:
		return PP_OUT 
	elif dbName == GU.CISBP:
		return CISBP_OUT
	else:
		return DEFAULT



####################################################################################################################
#This function will:
#1) Take the results of the motif comparison and use it to query the database of transcription factors
#2) Generate lists of refMotif object for each motif
#3) Add these to the masterObject for pickling
#4) Generate a file printout containing each motif and its TF matches (?)

#Where everything is updated and housed. Yea.
global masterMotifList
PP_DEFAULT_RETURN_KEYS = [GU.ID, GU.MOTIF, GU.TF, GU.SPECIES, GU.FAMILY, GU.DEATS]
##NOTE THAT FOR CISBP, the TF gene ID is in the "DBID" Column
CISBP_DEFAULT_RETURN_KEYS = ["DBID", "Motif_ID", "TF_ID", "TF_Species", GU.FAMILY, "DBDs"]  # Doesn't have a details column, has a DNA binding domains column instead

#TODO: HAVE THIS IMPLEMENT THE COMPARISON CALL AS WELL
if __name__ == '__main__':
	#Argument parsing
	parser = argparse.ArgumentParser(description='Examines the motif database to extract details about TF orthologs from HOMER comparison data..')
	parser.add_argument("comparisonPath", help="Path to the comparison data from HOMER compareMotifs.pl")
	parser.add_argument("-db", "--database", choices = ["PlantPan", "CISBP"], help="Which database to query, default is both.")
	parser.add_argument("-p", "--printOnly", action = "store_true", help = "Specify this if you only wish to generate an output file, not to update the master data list")
	parser.add_argument("-t", "--top", default = 5, help = "Up to how many orthologous binding site matches should we keep for each motif?  Default is 25- we want to keep most of them, many filtered out by species")
	parser.add_argument("-ms", "--matchScore", default=0.8, type=float, help = "Threshold match score for which motifs to keep.  Default is 0.8")
	parser.add_argument("-c", "--combined", action  = "store_true", help = "Specify if you wish to have a combined output file, default is separate")
	parser.add_argument("-s", "--species", help = "A list of species to search for in place of the default list  Format should be \"Genus species, Genus species...\"")  
	args = parser.parse_args()
	
	masterMotifList = pickle.load(open((os.getcwd() + "/v2" + GU.SER_MOTIF_MAP), "rb" ))
	#####Custom list of species to search
	if args.species:
		filterSpecies = [x.strip() for x in args.species.split(',')]
		if len(filterSpecies) == 0 or filterSpecies[0] == '' or filterSpecies[0] == ' ':
			print "Filter species specified are invalid.  Please try again."
			sys.exit()
	else:
		filterSpecies = GU.DEFAULT_SPECIES
	
	#Step through the files in the given directory:
	referenceMotifMap = dict() 
	if args.database: ##They specify which databse they wish to use
		if args.database == "PlantPan":
			referenceMotifMap = extractComparisonResults(args.comparisonPath + "/PP/homerResults", args.database, filterSpecies)
			
		elif args.database == "CISBP":
			referenceMotifMap = extractComparisonResults(args.comparisonPath + "/CIS-BP/homerResults", args.database,filterSpecies)
			
		writeRefMotifs(outFileName(args.database), referenceMotifMap)
	else: #Default, search both databses
		
		pp_map = extractComparisonResults(args.comparisonPath + "/PP/homerResults",GU.PP,  filterSpecies)
		cisbp_map = extractComparisonResults(args.comparisonPath + "/CIS-BP/homerResults", GU.CISBP,  filterSpecies)
		#If user wants the output files combined

		referenceMotifMap = combineReferenceMotifMaps(pp_map,cisbp_map)

		if args.combined:
			writeRefMotifs(outFileName("default"), referenceMotifMap)
			
		else: #Normal case, 2 separate files
			writeRefMotifs(outFileName(GU.PP),pp_map)
			writeRefMotifs(outFileName(GU.CISBP),cisbp_map)

	#Finally, update results and print results
	if not args.printOnly:
		updateMasterList(filterSpecies, masterMotifList, referenceMotifMap)
	
