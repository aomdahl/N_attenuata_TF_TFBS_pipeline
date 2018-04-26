#!/usr/bin/python
import sys
import os
import re
import argparse
import subprocess
from enum import Enum
sys.path.insert(0,"/home/likewise-open/ICE/aomdahl/Toolkit")
import GeneUtils as GU


#Map keys- from 
CS = "CandidateSequence"
DEATS = "Details"
EXACTS = "ExactMatchCount"
PVAL = "p-Value"
POS_FREQ = "PositiveFrequency"
NEG_FREQ= "NegativeBackgroundFrequency"
PMATRIX= "ProbabilityMatrix"
LODT = "LogOddsDetectionThreshold"
#True or false switch for if the motif appears in the overlap of orthologs
MOTIF_NOT_IN_OVERLAP = 0
MOTIF_IN_OVERLAP = 1

ORTHO_INDEX = 0
OVERLAP_STATE_INDEX  = 1




#This class stores each putative motif and all of the relevant, associated information
#Data members include:
	#Actual sequence
	#P-value
	#Genes in which it appears
	#Positive and negative frequencies with which it appears in the set
	#Conservation status [NoOrtholog | OrthologNoOverlap | OrthologAndOverlap] along with a conservation score.
	#Probability matrix
	#A list of all the genes containing the motif and their respective orthologs in Tomato, with a 0 or 1 indicating if the motif appears within that sequence
class pMotif:
	def __init__(self, mapData):
		self.motif = mapData[CS] #The actual sequence
		
		self.pvalue = mapData[DEATS][PVAL]
		self.posFreq = mapData[DEATS][POS_FREQ]
		self.negFreq = mapData[DEATS][NEG_FREQ]
		self.probMatrix= mapData[DEATS][PMATRIX]
		self.geneList = dict()  #GeneID: [OrthologID, 0/1 if Motif in their overlap]
		self.conservationType  = GU.NO_ORTHOLOG
		self.conservationScore = -1
		self.consPValue = -1
		#Either = or < or >, depending on what it is.  This is due to limited number of 
		#Bootstraps we can actually run, if its less than 1/1000
		self.consPValueCondition = ""
		
		self.logOdds  = mapData[LODT]#This is a value from HOMEr, we need to store it.
		
		#List of referenceMotifs that match this one
			#Not all pMotifs have this, so be sure to check the length
		self.refMotifs = list()
		self.blastPSearchDone = False
		
		
	def addGene(self, geneIn, orthoGene = "", motifInOverlap = MOTIF_NOT_IN_OVERLAP):
		self.geneList[geneIn] = [orthoGene, motifInOverlap]
		
	#To compare one pTFBS to another ("greater than")
	def __gt__(self, other):
		if self.pvalue > other.pvalue:
			if self.posFreq < self.posFre:
				return False
		if self.logOdds < other.logOdds:
			return False
		if self.conservationType != GU.ORTHO_AND_OVERLAP and other.conservationType == GU.ORTHO_AND_OVERLAP:
			return False
		if self.conservationScore < other.conservationScore:
			return False
		if self.consPValue > other.consPValue:
			return False
			
		return True
	
	def getGenes(self):
		return self.geneList
	
	def setConsPValue(self, value, condition):
		self.consPValue = value
		self.consPValueCondition = condition
		self.conservationType = GU.ORTHO_AND_OVERLAP
	
	#This tool returns a list of all genes associated with the motif that have orthologs.
	def getGenesWithOrthologs(self):
		retList = list()
		for gene in self.geneList:
			#If the gene has a valid ortholog already
			if self.geneList[gene][ORTHO_INDEX] != GU.EMPTY:
				retList.append(gene)
		return retList
		
	#Detailed annotations
	def addOrtholog(self, geneIn, orthologGene):
		
		#Check to make sure it hasn't been added already 
		if self.geneList[geneIn][ORTHO_INDEX] == GU.EMPTY:
			#print "adding the ortholog", orthologGene
			self.addGene(geneIn, orthoGene=orthologGene)
			
		if self.conservationType != GU.ORTHO_AND_OVERLAP:
				#And we assume that the gene has an ortholog if this call is being made at all.
				assert(self.geneList[geneIn][ORTHO_INDEX] is not None)
				self.conservationType = GU.ORTHO_NO_OVERLAP
	
	def setConservationType(self, c_type):
		self.conservationType = c_type
		
	#This function updates the overlap status for a gene in the current motif
	#Allows this to be done without knowing the ortholog explicitly, for convenience sake.
	#@param overlapStatus: 0 if not overlap, 1 if there is an overlap
	#@param gene: the gene in question
	
	def updateConservationExplicit(self, gene, overlapStatus):
		if overlapStatus == MOTIF_IN_OVERLAP:
			self.conservationType = GU.ORTHO_AND_OVERLAP
		else:
			#So the motif doesn't appear in the ortholog, 
			if self.conservationType != GU.ORTHO_AND_OVERLAP:
				#And we assume that the gene has an ortholog if this call is being made at all.
				assert(self.geneList[gene][ORTHO_INDEX] is not None)
				self.conservationType = GU.ORTHO_NO_OVERLAP
		
		self.geneList[gene][OVERLAP_STATE_INDEX] = overlapStatus
		
	#An alternative function, for use in other cases, to make sure that the pMotif's conservation status is up to date,
	#As best as is possible
	def updateConservation(self):
		if self.conservationScore != -1 or self.consPValue != -1:  #Do we have a conservation score?  If we do, conserved
			self.conservationType = GU.ORTHO_AND_OVERLAP #Do we have orthologs?  If so, we have orthologs
		elif len(self.geneList) >0:
			self.conservationType = GU.ORTHO_NO_OVERLAP
		else:  #Otherwise, no orthologs and no conservation.
			self.conservationType = GU.NO_ORTHOLOG	
		print "Updated to:", self.conservationType
	
	
	#Calculates the conservation score for each motif
	#This is done by:
	#    Number of conserved pairs of the motif
	# -----------------------------------------------
	# Num. of ortholog pairs where the  motif appears
	def calculateConservationScore(self):
		conservedPairs = 0
		orthologPairs = 0
		#Loop through list of genes,
		for gene in self.geneList:
			
			#If the gene has an ortholog...
			if self.geneList[gene][ORTHO_INDEX] == "":
				continue
			else:
				orthologPairs = orthologPairs + 1
			
			#If the motif overlaps between gene and ortholog
			if self.geneList[gene][OVERLAP_STATE_INDEX] == MOTIF_IN_OVERLAP:
				conservedPairs = conservedPairs + 1

		self.conservationScore = GU.calculateConservationScore(conservedPairs,orthologPairs)
		print "Calculated cons score", self.conservationScore
	
	#Allows you to get information about a given gene associated with a pMotif
	#@param geneID- the gene you are interested in
	#@param geneDataKey the info you want (i.e family, name, etc.)
	#@param refGene- True if you are interested in searching reference Genes, false if you want to search
	#@param actual ortholog genes.
	#@return None if no response,
	
	##TODO: Teset this out.  make sure it gives you what you want.
	def getGeneData(self, geneID, geneDataKey, refGene = True):
		if refGene:
			for ref in self.refMotifs:
				result = ref.findGeneData(geneID, geneDataKey)
				if result is not None:
					return result
			return None
		else:
			if self.refGeneNIATSites:
				geneID = geneID.replace("_", " ")
			if self.refGeneNIATSites[geneID]:
				return self.refGeneNIATSites[geneID]
		return None
				
###########################BlastP Lookup for member reference Motifs:
	#Calls the Blastp search on all the member reference motifs, and stores the results in a hand data structure:
	#  NIATGeneID : [Species1GeneSet, Species2GeneSet, Species3GeneSet...]
	#@param consPValue: conservation p-value threshold for this
	#@param return count, how many hits to accept on each blast search
	#@param eThresh, the e-value threshold
	#@consPVal: conservation pValue to accept
	def BLASTpSearchRefMotifs(self, returnCount =3, eThresh =-1):
		self.blastPSearchDone = False
		#this is the storage lcoation for all the reference genes, organized by NIAT gene hit
		#Subkeys include
		self.refGeneNIATSites = dict()
		
		#List keeping track of gene search hits and their associated e-values, for write out.
		self.blastEvals = dict()
		self.blastPercentIdentity = dict()
		
		print "Performing the BLASTp search for genes enriched with the motif", self.motif
		
		for rm in self.refMotifs:
			hitList = rm.BLASTPSearchRefs(returnCount, eThresh)
			self.addRefGeneNIATSites(hitList, self.refGeneNIATSites, self.blastEvals, self.blastPercentIdentity)
		self.blastPSearchDone = True
		return self.refGeneNIATSites
		
	#This tool actually updates the refGeneNIATSites object
	#This object is a map of NIAT gene IDs -->Species -->SpeciesGenes
	#@param referenceList: the results of the BlastP search, a table of [speciecs, query, niat_hit, %identity, evalue] lists
	#@param refMap to update
	#@param evalueMap: list of evalues associated with the reference genes
	#@param percIDMap: list of percentage blast identity  associated with the reference genes
	def addRefGeneNIATSites(self, refList, refMap, evalueMap, percIDMap):
		
		if not refList or len(refList) == 0:
			return
	
		for match in refList:
			species = match[0]
			query = match[1]
			niat_hit=match[2]
			thresh=match[3]
			e_val = match[4]

			if niat_hit not in refMap:
				refMap[niat_hit] = dict()
			if species not in refMap[niat_hit]:
				refMap[niat_hit][species] = set()
			#The actual adding call
			#print "Adding into our repository", query, niat_hit
			refMap[niat_hit][species].add(query)
			
			#Update the eval as well, for write out.
			evalueMap[niat_hit] = e_val
			#Update the percentage overlap.  Need this for output writing.
			percIDMap[niat_hit] = thresh

			
###############################Output stuff, mostly for debugging #########################

	#File write out for a given motif.
	#Modify this.
	def writeOut(self, outFile):
		outFile.write(">"+self.motif + '\t' + self.pvalue + '\t' + self.posFreq + '\t' + self.negFreq)
		if self.conservationType == GU.ORTHO_AND_OVERLAP or self.conservationType == GU.ORTHO_NO_OVERLAP:
			outFile.write("\nConservation Status: " +self.conservationType + '\n')
			if self.conservationType == GU.ORTHO_AND_OVERLAP :
				outFile.write("Conservation Score = " + str(self.conservationScore) + '\t')
				outFile.write("Conservation custom p-value " + self.consPValueCondition + " " + str(self.consPValue))
			outFile.write("\n")
			self.writeGenesOrderedByConservation(outFile)
		outFile.write(self.probMatrix+ '\n')
	
	#This will write out a file in Homer format.
	#@param outfile: the file to write out to.
	def writeOutHomer(self, outFile):
		outFile.write(">"+self.motif + '\t' + str(self.pvalue) + ":" + self.motif + '\t' + str(self.logOdds) + '\t' + str(self.conservationScore) + '\t' + str(self.consPValue) + '\n')
		outFile.write(self.probMatrix)
	##
	
	#This function aids in printing out the list of genes in a meaningful order
	#Prints in the order of those with overlaps, those with orthologs, genes with no orthologs.
	#@param outfile
	def writeGenesOrderedByConservation(self, outFile):
		#Iterate through the map 3 times, each time under different conditions
		#First: print all those with Overlap 
		for gene in self.geneList:
			if self.geneList[gene][OVERLAP_STATE_INDEX] == MOTIF_IN_OVERLAP:
				outFile.write(gene + ":" + self.geneList[gene][ORTHO_INDEX]+ " Overlap Status:" + str(self.geneList[gene][OVERLAP_STATE_INDEX]) + '\n')
	
		#Second: print all those with orthologs, no overlap
		for gene in self.geneList:
			if self.geneList[gene][OVERLAP_STATE_INDEX] == MOTIF_NOT_IN_OVERLAP and self.geneList[gene][ORTHO_INDEX] != GU.EMPTY:
				outFile.write(gene + ":" + self.geneList[gene][ORTHO_INDEX]+ " Overlap Status:" + str(self.geneList[gene][OVERLAP_STATE_INDEX]) + '\n')
		#Finally: print all those with no 
		for gene in self.geneList:
			if self.geneList[gene][ORTHO_INDEX] == GU.EMPTY and  self.geneList[gene][OVERLAP_STATE_INDEX] == MOTIF_NOT_IN_OVERLAP:
				outFile.write(gene + ":" + self.geneList[gene][ORTHO_INDEX]+ " Overlap Status:" + str(self.geneList[gene][OVERLAP_STATE_INDEX]) + '\n')
