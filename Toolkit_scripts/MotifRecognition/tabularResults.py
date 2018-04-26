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
import glob
from MotifList import MotifList
from refMotif import refMotif
import errorReporting as ER
sys.path.insert(0,"/home/likewise-open/ICE/aomdahl/Toolkit")
import HomerTools as HT
import GeneUtils as GU
from getModuleOrthologs import parseNIATMotifLocations
from getModuleOrthologs import geneOrtholog
	
####################Write out conservation results




####################Write out TF results##############################

def writeTFResultsTable(output, referenceSpecies):
	outfile = open(output, 'w')
	
	########## Write the header out ########
	outfile.write("Motif\tMotifMatchScore\tpValue\tFrequencyInSet\tReferenceMotif\t")
	#Write the species headers
	referenceSpecies = masterMotifList.referenceSpecies
	for species in referenceSpecies:
		outfile.write(species.replace(" ", "_"))
		outfile.write('\t')
	outfile.write('\n')
	######### Write the sections

	writeResultsPreBLAST(outfile, masterMotifList.getConserved(), referenceSpecies, unique = True)
			
#Query a hash table to see if contains a gene
#if it doesn't we add it to it and return true
#else, return false.
def uniquenessChecker(hashTable, geneID):
	if geneID == "NA":
		return True
	if geneID in hashTable:
		return False
	else:
		hashTable[geneID] = 1
		return True

#This method will indicate if a header should be written or not on a file
#i.e. no if the file exists already, and then writes it
#@param appendStatus: true if we appendign to files, otherwise just overwrites
#@param filePath to possibly write to
#@param headerList to write out
def headerAppendWrite(appendStatus, filePath, headerList):
	if appendStatus:
		return
	else:
		tempOpen = open(filePath, 'w')
		GU.writeTabList(tempOpen, headerList)
		tempOpen.close()

##Prints a report mapping motifs to identified TF in other species.
#@param outputFile- fileStream to write out to
#@param pTFBS that we are currently editing
#@param unique: only show unique genes-i.e. show a gene only once.
def writeResultsPreBLAST(outputFile, pTFBSList,referenceSpecies, verbosity = 1, unique = True):
	for pTFBS in pTFBSList:
		if pTFBS.consPValue == -1 or  pTFBS.consPValue > args.consPval:
			continue  #onto the next iteration, we aren't interested in writing out about this pTFBS
		writeOutTable = list()
		#Get the list of related TF genes
		refMotifs = pTFBS.refMotifs 
		uniqueGenes = dict() #Keep track of genes that have already been reported so we only report them once.

		if len(refMotifs) <= 0: #Case with no reference motifs, need a blank line with the placeholders
			writeTable = [pTFBS.motif,"-", pTFBS.pvalue,pTFBS.posFreq, ""]
			outputFile.write("\t".join(writeTable) + '\n')
			continue
		#Loop through the list of referenceMotifs/TFs
		for referenceMotif in refMotifs:
			speciesRefList = list()
			matchScore = referenceMotif.findGeneData(referenceMotif.name, GU.MATCH_SCORE)
	
			#Get a list for each species
			for species in referenceSpecies:
				speciesRefList.append(referenceMotif.getFromSpecies(species))
			
			#print "Current focus hit is"
			#raw_input(referenceMotif.name)
			#Loop through each list until all items have been visited
			longestListLength = max(len(i) for i in speciesRefList)
			#print "Length of longest list"
			#raw_input(longestListLength)
			#referenceMotif.printOut()
			#raw_input("...")
			for i in range(0, longestListLength):
				writeTable = [pTFBS.motif,matchScore, pTFBS.pvalue,pTFBS.posFreq,referenceMotif.name]
				dataExtendList = list()
				#Loop through each sublist
				for subList in speciesRefList:
					if i < len(subList):
						if unique:
							if uniquenessChecker(uniqueGenes, subList[i][GU.ID]): #If we have a unique item, append it to the list
								appendGeneData(subList[i], dataExtendList)
							else: #We don't 
								#Can we grab the next option?-- work on this, but at a later time.
								possibleAlternative = False
								for j in range(i, len(subList)):
									if uniquenessChecker(uniqueGenes, subList[j][GU.ID]):
										appendGeneData(subList[j], dataExtendList)
										possibleAlternative = True
										break;
								if not possibleAlternative:
									dataExtendList.append("-")
						else:
							appendGeneData(subList[i], dataExtendList)
					else:
						dataExtendList.append("-")
						
				if all(x == "-" for x in dataExtendList):
					#print "Empty line", dataExtendList
					pass
				else:
					finalList = writeTable + dataExtendList
					#print "Final list",  finalList
					outputFile.write("\t".join(finalList)+ '\n')
				
				writeTable = [pTFBS.motif,matchScore, pTFBS.pvalue,pTFBS.posFreq,referenceMotif.name]
				
#Returns a table with additional headers for write out, based on the verbosity
def verbosityMetrics(pTFBS, verbosity):
	if verbosity == 1:
		return [pTFBS.motif, str(pTFBS.pvalue), str(pTFBS.posFreq)]
	if verbosity == 2:
		return [pTFBS.motif, str(pTFBS.pvalue),  str(pTFBS.posFreq), str(pTFBS.consPValue)] 
	if verbosity >= 3:
		return [pTFBS.motif, str(pTFBS.pvalue),  str(pTFBS.posFreq), str(pTFBS.negFreq), str(pTFBS.conservationScore), str(pTFBS.consPValue)] 

	

#Returns a string for the header based on the verbosity
def verbosityHeader(speciesList, verbosity):
	headerString = ""
	if verbosity is 1:
		headerString = "Motif\tpValue\tFrequencyInSet\tPutativeNIATTFGene\t"
	if verbosity is 2:
		headerString = "Motif\tpValue\tFrequencyInSet\tConservationPVal\tPutativeNIATTFGene\t"
	if verbosity is 3:
		headerString = "Motif\tpValue\tFrequencyInSet\tNegativeFrequency\tConservationScore\tConservationPVal\tPutativeNIATTFGene\t"

	#Write the species headers
	for s in speciesList:
		headerString = headerString + (s.replace(" ", "_")) + '\t'
		if verbosity >2:
			headerString = headerString + "%Identity\t"
	return headerString + '\n'


#This uses the internal components of the MotifList structure to print out the results of the blastp saerch
#THIS CAN ONLY BE CALLED AFTER A BLAST PSEARCH.
#@param outputFile to write to
#@param pTFBSList
#@param verbosity- how much to say

#Bioinformaticicst results:
#Module gene ID
#Module name
#Module motifs
#Motif conservation score, appearance frequency, etc. etc.
#Essentially a modification of what you ahve above actually.
def writeResultsPostBlast(outputFile, pTFBSList, speciesList, verbosity):
	outfile = open(outputFile, 'w')
	outfile.write(verbosityHeader(speciesList, verbosity))
	
	for pTFBS in pTFBSList:
		#Only keep those with a statistically significant conservation score.
		if pTFBS.consPValue == -1 or  pTFBS.consPValue > args.consPval:
			continue
			
		writeOutTable=list()
		writeOutTable = verbosityMetrics(pTFBS, verbosity)
		#Make sure the search has actually been made.  If not, we have a houston problem!
		try:
			geneHitList = pTFBS.refGeneNIATSites
			refMap[niat_hit][species].add(query)
		except Exception:  #In the event that a certain NIAT hit isn't actually in the reference map, do the search
			try:
				if not pTFBS.blastPSearchDone: #Also may be the case that the search was never done, may not have the variable assigned.
					pTFBS.BLASTpSearchRefMotifs()
					geneHitList = pTFBS.refGeneNIATSites
			except AttributeError:
				pTFBS.BLASTpSearchRefMotifs()
				geneHitList = pTFBS.refGeneNIATSites
				
		if len(geneHitList) == 0:
			print "No gene hits for ", pTFBS.motif, " motifs"
			print "Please validate this result."
			for s in speciesList:
				writeOutTable.append('-')
			outfile.write("\t".join(writeOutTable))
			outfile.write('\n')

		for n_gene in geneHitList:
			#Create the list keeping track of list iteration
			print "This gene is in the hit list...", n_gene
			iterList = list()

			for species in speciesList:
				if species in geneHitList[n_gene]:
					iterList.append(list(geneHitList[n_gene][species]))
				else:
					iterList.append([]) #If the species has no hits, then just append an empty list.

			longestListLength = max(len(i) for i in iterList)

			for i in range(0, longestListLength):
				writeOutTable.append(n_gene)
				for s in iterList:
					if i < len(s):
						familyGroup = pTFBS.getGeneData(s[i], GU.FAMILY)
						if familyGroup:
							writeOutTable.append(s[i] + "," + familyGroup)
						else:
							writeOutTable.append(s[i] +  ",NA")
					else:
						writeOutTable.append("-")

				outfile.write("\t".join(writeOutTable))
				outfile.write('\n')
				writeOutTable = verbosityMetrics(pTFBS, verbosity)
	outfile.close()
		

#This helper funtion adds the gene data to our matrix
#Pulled out into a function to help make the code a bit cleaner
#@param refMotifMap: map from a reference motif for a certain gene TF homolog
#@param writeTable: table to append to for writing out.
def appendGeneData(refMotifMap, writeTable):
	if refMotifMap[GU.FAMILY] == " " or refMotifMap[GU.FAMILY] == "":
		writeTable.append(refMotifMap[GU.ID] + ',' + refMotifMap[GU.TF])
	else:
		writeTable.append(refMotifMap[GU.ID] + ',' + refMotifMap[GU.FAMILY])
	

#Module-based results: this outputs results for each gene in the module, where results exist
#ouput table contains the following:#
#1)List all of the genes in the current module
#2)List their functional annotation (gff file lookup)
#3)List which module motifs appear in the promoter region of these genes####WROGN
#4)Give the conservation score of that motif throughout the module
#5) List the conservation pvalue of the motif

def moduleResultsLibrary(masterMotifList, append = True, outpath = "Module_lib.tsv"):
	writeState = 'w'
	if append:
		writeState = 'a'
		
	fileStream = open(outpath, writeState)
	#Get the list of genes in the motif
	modGeneList = GU.getModuleGenesList(os.getcwd())
	#Get the motifs in the promoters of the different genes
	geneToMotifMap = loadGenePromoters()
	if not modGeneList:
		modGeneList = GU.getModuleGenesList(os.pardir)
		if not modGeneList:
			print "Error: no module genes file available.  Please ensure a module file is in the local directory."
			return
	
	for gene in modGeneList:
		annot = GU.getFunctionalAnnotations(gene)
		try:
			moduleName = masterMotifList.moduleName
		except AttributeError:
			moduleName = "None"
	#moduleName = masterMotifList.moduleName  for now, just override

		try:
			motifList = geneToMotifMap[gene].motifList
		except KeyError:
			#"No module motifs appear in this gene's promoter sequence"
			motifList = ["NA"]
			GU.writeTabList(fileStream, [gene, annot, "NA", moduleName, "NA", "NA"])
			continue
		
		for promoterMotif in motifList:
			consScore = masterMotifList.getMotif(promoterMotif).conservationScore
			consPVal = masterMotifList.getMotif(promoterMotif).consPValue
			GU.writeTabList(fileStream, [gene, annot, promoterMotif, moduleName, consScore, consPVal])
			
	fileStream.close()

	
	
	
#Helper function for moduleResultsLibrary
#this makes the call to identify the motifs found in the promoter sequences of the differnet genes.
def loadGenePromoters():
	if not os.path.isfile("topMotifSites.bed"):
		print "It appears that the necessary gene : motif BED file mapping isn't available"
		print "Please run scanMotifGenomeWide.pl"
		print "i.e. scanMotifGenomeWide.pl topMatches.motif 1kb_promoters.fa -bed > topMotifSites.bed"
		return
	else:
		print "Opening topMotifSites.bed file"
		print "Module Results Library written out"
		bedResults = parseNIATMotifLocations("topMotifSites.bed", updateSuper = False)
		return bedResults[0]

#This creats a .tsv file summarizing potential TFs correlated with the module genes.  Headers include the following data
#1)Identified potential TF gene
#2) e-value for that gene on the blast search
#3) The gene's functional annotation, if available
#2)pTFBS for that gene
#3)The match score of the predicted matrix with the given motif
#3)The module it is predicted to regulate
#3)The gene's correlation with that module (correlation of its expression to average expression of module genes)
#4) The module with which it is coexpressed
#5) Its connectivity with that module
#@param moduleName- the module of the current putativeTranscription factor
#@param motifObject- the motif Object associated with the motif
#@param fileStream to write out to
#Potential spot to speed up by making that call only once.
def pTFLibrary(moduleName, motifObject, fileStream,showMissing = False):
	
	print "Writing a putative TF library..."
	#Get a list of correlations for the gene with the module
	geneToModuleCorrelation = GU.expressionCorrelationLookup(motifObject.refGeneNIATSites.keys(), GU.getModuleGenesList("./", pathOnly = True))
	print "Length of potential alternatives:", len(motifObject.refGeneNIATSites)
	for ref in motifObject.refGeneNIATSites:
		if ref == "No hits" and showMissing:  #Case where have not mapped BLASTpgenes
			continue
		#Ref is the NIAT gene
		funct = GU.getFunctionalAnnotations(ref)
		bindingMotif = motifObject.motif
		regulatingModule = moduleName
		expressionModule = GU.getGeneModule(ref, scaledConnectivity = True)
		matchScore = motifObject.getGeneData(ref, GU.MATCH_SCORE, True)
		try:
			evalue = motifObject.blastEvals[ref]
			percentIdentity = motifObject.blastPercentIdentity[ref]
		except KeyError:
			evalue = "NA"
			percentIdentity = "NA"
		try:
			moduleCorrelation = geneToModuleCorrelation[GU.removeGeneSuffix(ref)]
		except KeyError:
			moduleCorrelation = "Data unavailable"
			print "Unable to determine modular correlation for:", ref
			
			#If verbosity is high, include all data
		if args.verbosity > 1:
			GU.writeTabList(fileStream, [ref,evalue, funct, percentIdentity, bindingMotif, regulatingModule, moduleCorrelation, expressionModule[0], expressionModule[1]])
		elif ref in GU.FALSE_RESULTS or ref == "":
			continue
		else:  #If its only 1, only record the most significant ones.
			print "Module for gene", ref, "has a correlation of", moduleCorrelation, "with the current gene module"
			if moduleCorrelation  in GU.FALSE_RESULTS:
				if ref not in GU.FALSE_RESULTS:
					GU.writeTabList(fileStream, [ref,evalue, percentIdentity, funct, bindingMotif, matchScore, regulatingModule, moduleCorrelation, expressionModule[0], expressionModule[1]])
				continue
			if float(moduleCorrelation) >= args.corr or float(moduleCorrelation) <= (args.corr * -1):
				#print "Module for gene", ref, "has a correlation of", moduleCorrelation, "with the current gene module"
				GU.writeTabList(fileStream, [ref,evalue, percentIdentity, funct, bindingMotif,matchScore, regulatingModule, moduleCorrelation, expressionModule[0], expressionModule[1]])

#This prints out a report on all of the motifs found
#Specifically, their scores, conservation, etc. etc.
#@param masterMotifList- the big list of motifs
#@param outpath -file path to outpuut
def motifReport(masterMotifList, outpath):
	
	outStream = open(outpath, 'a')
	for pTFBS in masterMotifList.motifMap:
		cm = masterMotifList.motifMap[pTFBS] #cm stands for current motif
		GU.writeTabList(outStream, [cm.motif, cm.pvalue, cm.posFreq, cm.negFreq, cm.conservationType, cm.conservationScore, cm.consPValue])
	outStream.close()

			
	
			
#Wrapper function for library printing tools
#Trying to be a good coder and not repeat code.
#@param masterMotifList the serialized list
def libraryPrint(masterMotifList, libMethod, outpath, append = False):
	try:
		moduleName = masterMotifList.moduleName
	except AttributeError:
		moduleName = "None"
		
	appendStatus = 'a'
	if not append:
		appendStatus = 'w'
		
	if libMethod == "pTF":
		
		headerAppendWrite(append, outpath, ["PutativeTF","e-Value", "BLASTPercentIdentity", "FunctionalAnnotation", "Put.BindingMotif", "MotifLibraryMatchScore","RegulatingModule", "ModuleCorrelation", "ExpressionModule", "ModuleExpressionScore"])
		outFile = open(outpath, 'a')
		for pTFBS in masterMotifList:
			motifObj = masterMotifList.getMotif(pTFBS)
			if motifObj.consPValue <= args.consPval and motifObj.consPValue != -1:
				pTFLibrary(moduleName, motifObj, outFile)
			
		print "Putative TF library written to", outpath
		
	elif libMethod == "Motifs":
		headerAppendWrite(append, outpath, ["Motif","p-Value", "%FreqInPromoters", "%FreqInBackground", "ConservationType","ConservationScore", "ConservationP-value"])
		motifReport(masterMotifList, outpath)
		print "Motif record written out to", outpath
		
	elif libMethod == "moduleResultsLib":
		headerAppendWrite(append, outpath,["ModuleGene", "FunctionalAnnotation", "PromoterBindingMotif", "ModuleName", "MotifConservationScore", "p-value"])
		moduleResultsLibrary(masterMotifList, append = True, outpath=outpath)
		

	


####################################################################################################################
#This script:
#1) Access the giant dump file we have created
#2) Parse its data and print it out into a customizeable table for human viewing.
global masterMotifList


if __name__ == '__main__':
	#Argument parsing
	parser = argparse.ArgumentParser(description='Examines the motif database to extract details about TF orthologs from HOMER comparison data..')
	parser.add_argument("-o", "--output", help="Output File Name and path")
	parser.add_argument("-v", "--verbosity", choices = [1,2,3], type = int, default =1, help="Specify how much detail level to include in output reports.  Default is 3, the most information" + \
						"Note that selecting 1 will only output correlating gene families")
	parser.add_argument("-s", "--species", help = "The species to include in the report: \"Genus species, Genus species\".  Default is standard list")
	parser.add_argument("-t", "--reportType", choices=["PostBlast", "RefMotifs", "pTFLib", "moduleResultsLib", "PreBlast", "Motifs"], default = "RefMotifs", help ="Please specify what kind of report you want: a post-blast report, a report of reference motifs, etc")
 	parser.add_argument("-c", "--conservations", choices = ["all", "orthologs", "none", "conserved"], help = "Choose this option to specify and generate the gene lists that were either orthologs, conserved etc.")
 	parser.add_argument("-cp", "--consPval", type = float, default = 0.01, help = "Specify the threshold conservation p-value to retain in this report.")
	parser.add_argument("--corr", type = float, default = 0.5, help = "Specify the intramodular correlation of putative TF genes within the genes they putatively regulate, default is 0.5")
	parser.add_argument("-a", "--append", action="store_true", help = "Select this option if you wish to append to files that already exist.  Default is false")
	args = parser.parse_args()
	
	#TODO: condense this list- its reduntant with the library function call above.
	speciesList = GU.DEFAULT_SPECIES
	if args.species:
		speciesList = args.species
	
	try:
		if args.conservations:
			writeOutConservation(args.conservations, args.verbosity)
		
		if args.reportType != "PreBlast":
			masterMotifList = pickle.load(open(os.getcwd() + "/vF" + GU.SER_MOTIF_MAP, "rb" ))

		if args.reportType == "PostBlast":
			writeResultsPostBlast(args.output, masterMotifList.getConserved(), speciesList, args.verbosity)
			print "Wrote out post-blast search results to", args.output
		
		elif args.reportType == "pTFLib":
			if args.verbosity >1:
				libraryPrint(masterMotifList, "pTF", args.output, append = args.append) ##??

			else:
				libraryPrint(masterMotifList, "pTF", args.output, append = args.append)
		#In the case of a motifs report....
		elif args.reportType == "Motifs":
			output = args.output
			if not args.output:
				output = "MotifHits.tsv"
			libraryPrint(masterMotifList, args.reportType, output)
		
		elif args.reportType == "moduleResultsLib":
			output = args.output
			if not args.output:
				output = "Module_lib.tsv"
			libraryPrint(masterMotifList, args.reportType, output)

		elif args.reportType == "PreBlast":
			masterMotifList = pickle.load(open(os.getcwd() + "/v3" + GU.SER_MOTIF_MAP, "rb" ))
			writeTFResultsTable(args.output, speciesList)
			print "Wrote out pre-blast search results to", args.output
		
		else:
			print "Please specify a report type"
			sys.exit()
	except IOError:
		print "It appears that the necessary serialized data file is missing."
		print "Please ensure you've specified the correction options, are in the correct directory and have run previous analysis already"
		sys.exit()
