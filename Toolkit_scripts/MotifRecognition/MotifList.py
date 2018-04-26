#!/usr/bin/python
import sys
import os
import re
import argparse
import subprocess
from enum import Enum
from pMotif import pMotif
sys.path.insert(0,"/home/likewise-open/ICE/aomdahl/Toolkit")
import GeneUtils as GU
from refMotif import refMotif

#Map keys- from 
CS = "CandidateSequence"
DEATS = "Details"
EXACTS = "ExactMatchCount"
PVAL = "p-Value"
POS_FREQ = "PositiveFrequency"
NEG_FREQ= "NegativeBackgroundFrequency"
PMATRIX= "ProbabilityMatrix"

#This class stores the list of Motifs for easy data access.
#
class MotifList:
	
	def __init__(self):
		self.motifMap = dict() #Key is the sequence, value is the pMotif object
		self.iterLocation = 0
		#Used for iterator.
		self.currentMotif = None
		self.referenceSpecies = list()
		self.moduleName = ""
	
	def addMotif(self, motifObj):
		self.motifMap[motifObj.motif] = motifObj

	def getMotif(self, motifSeq):
		return self.motifMap[motifSeq]

#Tool to merge two lists together if needed
#@param list to merge
	def append(self, motifListObj):
		#Merge the pTFBS list
		for pTFBS in motifListObj.motifMap:
			if pTFBS in self.motifMap:
				print "The motif", pTFBS, "is already in this map.  The one with the highest signficance and conservation will be kept"
				if self.motifMap[pTFBS] > motifListObj.motifMap[pTFBS]:
					pass
				else: #Reassig the newer one.
					self.motifMap[pTFBS] = motifListObj.motifMap[pTFBS]  
								
			else:
				self.motifMap[pTFBS] = motifListObj.motifMap[pTFBS]
				
		#update the motifList as well., only the unique items
		for element in motifListObj.referenceSpecies:
			if element not in self.referenceSpecies:
				self.referenceSpecies.append(element)
		
		
#Not sure if this will work, can try out someother time.
	def __iter__(self):
		return self
	
	def next(self):
		if self.iterLocation < len(self.motifMap):
			self.currentMotif = self.motifMap.keys()[self.iterLocation]
			self.iterLocation += 1
			return self.currentMotif
		else:
			self.iterLocation = 0
			raise StopIteration()
			
			
	#Calculate the conservation scores for each gene in the list:
	def calculateConservationScores(self):
		for motifSeq in self.motifMap:
			self.motifMap[motifSeq].calculateConservationScore() 

	
	#This is a shortcut method for adding an ortholog to the list of motifs and their associated genes.
	#@param origina_gene: the gene we are adding an ortholog to
	#@param ortholog: the ortholog for that gene
	def addOrtholog(self, original_gene, ortholog):
		for motifSeq in self.motifMap:
			if original_gene in self.motifMap[motifSeq].geneList:
				self.motifMap[motifSeq].addOrtholog(original_gene, ortholog) #This automatically updates the conservation type of the motif as well.
	
	#This implicitly updates the conservation status for each pTFBS based on its data members
	def updateConservationStatuses(self):
		for motifSeq in self.motifMap:
			self.motifMap[motifSeq].updateConservation()
			
	###########These methods sort the motifs by their conservation status###################

	#Returns a list of all motifs that are conserved, sorted by their conservation score.
	def getConserved(self):
		#print self.motifMap
		#self.updateConservationStatuses()
		retMotifs = list()
		for motifObject in self.motifMap:
			if self.motifMap[motifObject].conservationType == GU.ORTHO_AND_OVERLAP:
				retMotifs.append(self.motifMap[motifObject])
			else:
				print self.motifMap[motifObject].conservationType
		return sorted(retMotifs, key = lambda element:element.conservationScore)
	
	def getSigConserved(self, conservationP = 0.001):
		self.updateConservationStatuses()
		retMotifs = list()
		for motif in self.getConserved():
			if motif.consPValue != -1 and motif.consPValue <= conservationP:
				retMotifs.append(motif)
		return retMotifs
		
	#Returns a list of all motifs with no orthologs at all
	def getNoOrtholog(self):
		self.updateConservationStatuses()
		retMotifs = list()
		for motifObject in self.motifMap:
			if self.motifMap[motifObject].conservationType == GU.NO_ORTHOLOG:
				retMotifs.append(self.motifMap[motifObject])
		return retMotifs
	
	#Returns a list of all motifs with orthologs but no conservation
	def getOrthologNoOverlap(self):
		retMotifs = list()
		for motifObject in self.motifMap:
			if self.motifMap[motifObject].conservationType == GU.ORTHO_NO_OVERLAP:
				retMotifs.append(motifObject)
		return retMotifs
		
	
	###################Call all the pTFBS#############
	
	#Do a blast search on each member of the list
	#Store and update each PTFBS in the process.
	def BLASTpSearchAllMotifs(self, returnCount , eThresh, consPValue = 1):
		self.blastPSearch = False
		print "Performing a BLASTp search on all putative motifs..."
		for motif in self.motifMap:
			#Do we want to limit this search by pvalue
			if float(self.motifMap[motif].pvalue) > consPValue or self.motifMap[motif].pvalue == -1:  #-1 if the results not initialized
				print "Not doing for some reason..."
				print "Must be smaller than..", consPValue
				print "pval", self.motifMap[motif].pvalue
 				continue
			else:
				self.motifMap[motif].BLASTpSearchRefMotifs(returnCount, eThresh)
		print "BLASTp search complete."
		self.blastPSearch = True

	
	#########Tools for writing to a file the data in the lists###################################
	
	#Writes all three lists
	#@param: paths for all three categories
	def writeConservationLists(self, pathConserved, pathNoOrtholog, pathNoOverlap):
		writeConserved(pathConserved)
		writeNoOrtholog(pathNoOrtholog)
		writeOrthologNoOverlap(pathNoOverlap)

	def writeConserved(path):
		writeList(path, self.getConserved())
	
	def writeNoOrtholog(path):
		writeList(path, self.getNoOrtholog())
	
	def writeOrthologNoOverlap(path):
		writeListHomerFormat(path, self.getOrthologNoOverlap())
	
	def writeSigConservedHomer(path):
		writeList(path, self.getSigConserved())
		
	def writeList(self, path, c_list):
		motifFile = open(path, 'w')
		for motif in c_list:
			motif.writeOut(motifFile)
		motifFile.close()
				
	def writeListHomerFormat(self, path, c_list):
		motifFile = open(path, 'w')
		for motif in c_list:
			motif.writeOutHomer(motifFile)
		motifFile.close()

	
	#Write the list out to a file for checking purposes.
	def writeToFile(self, filePath):
		outFile = open(filePath, 'w')
		for motifObj in self.motifMap:
			self.motifMap[motifObj].writeOut(outFile)
		outFile.close
		

