#!/usr/bin/python
import sys
import os
import re
import argparse
import subprocess
from enum import Enum
import GeneUtils as GU
sys.path.insert(0,"/home/likewise-open/ICE/aomdahl/Toolkit")
import HomerTools as HT
import BLASTPSearch #import searchBLASTPDatabase



#This class stores information regarding reference motifs/TFS and are ALWAYS TO BE ASSOCIATED# with a pMotif object.
class refMotif(object):
	#@param db_entry_maps: a list containing maps of TF matches.  
	##These maps are characterized by the following:
	
		#motif = db_entry_map[GU.MOTIF]
		#TFName = db_entry_map[GU.TF]
		#species = db_entry_map[GU.SPECIES]
		#family = db_entry_map[GU.FAMILY]
		#details = db_entry_map[GU.DEATS]
	
	def __init__(self, db_entry_maps, db):#Hold up. 
		self.name = db_entry_maps[0][GU.MOTIF] #The name really should be the same for all of them.... hopefully
		self.mapList = db_entry_maps
		
		for entry in self.mapList:
			entry[GU.SPECIES] = entry[GU.SPECIES].replace("_", " ")
		#consider special cases for CISBP database.
		self.db = db
		
	def setMatchParams(self, rank, score):
		self.matchScore = score
		self.matchRank = rank
		
	
	#Write the refMotif out to a file
	#@param openFile stream
	def writeToFile(self, openFile):
		openFile.write(self.name + " Match Score: " + self.matchScore + '\n')
		for r in self.mapList:
			openFile.write(r[GU.ID]+ "\t" + r[GU.TF]+ "\t" + (r[GU.SPECIES]).replace("_", " ") + '\t' + r[GU.FAMILY])
			if len(r[GU.DEATS]) > 20:
				openFile.write('\n' + r[GU.DEATS] + '\n')
			elif len(r[GU.DEATS]) > 0:
				openFile.write('\t' + r[GU.DEATS] + '\n')
			else:
				openFile.write('\n')
			
		openFile.write('\n')
		
	#Print it to the console:
	#@param openFile stream
	def printOut(self):
		print self.name
		#raw_input(self.name)
		for r in self.mapList:
			print(r[GU.ID]+ "\t" + r[GU.TF]+ "\t" + (r[GU.SPECIES]).replace("_", " ") + '\t' + r[GU.FAMILY])
			if len(r[GU.DEATS]) > 20:
				print('\n' + r[GU.DEATS])
			elif len(r[GU.DEATS]) > 0:
				print ('\t' + r[GU.DEATS])
			else:
				print

	
	#A simple helper function that identifies the data for  a given searched gene
	#@param the reference gene to search for
	#@param the key to search for (i.e. family, name, motif, etc. etc.)
	#@return result or None if there is none.
	def findGeneData(self, refID, key):
		if key == GU.MATCH_SCORE:
			return self.matchScore
		if key == GU.MATCH_RANK:
			return self.matchRank 
		for subRef in self.mapList:
			if (refID in subRef[GU.ID]) or (subRef[GU.ID] in refID) :
				return subRef[key]
		
		return None
	
	#Returns a list of all of the TF matches for this particular motif for a specified species
	#@param speciesName to search for
	#@return a list of all the maps for a given species
	def getFromSpecies(self, speciesName):
		retList = list()
		#print "Name of current match motif"
		#raw_input(self.name)
		for tf in self.mapList:
			if tf[GU.SPECIES] == speciesName:
				retList.append(tf)
			elif tf[GU.SPECIES].replace("_", " ") == speciesName:
				retList.append(tf)
			else:
				continue
		return retList
	
	#Helper function that returns a list of the species in this particular reference motif
	def getSpecies(self):
		speciesList = list()
		for tf in self.mapList:
			if tf[GU.SPECIES] not in speciesList:
				speciesList.append(tf[GU.SPECIES])
		return speciesList
	
	#Gets lists of TF refs sorted by species
	#@return a dictionary of species: the gene info.
	def sortedBySpecies(self, genesOnly = False):
		masterList = dict()
		speciesList = self.getSpecies()
		for species in speciesList:
				if genesOnly:
					masterList[species] = [i[GU.ID] for i in self.mapList if (i[GU.SPECIES] == species)]
				else:
					masterList[species] = [i for i in self.mapList if (i[GU.SPECIES] == species)]
		return masterList  #Keep those with the correct query search

			
	
	#This performs the blastP search on all of the elements of the refMotif object
	#@param returnCount- how many hits to take for each one
	#@param eThresh: the e-value cutoff.
	#@return a list of all the hits in the following format
	#[ [GeneID, GeneSpecies, %identity, e-value], [GeneID, GeneSpecies...]...]
	def BLASTPSearchRefs(self, returnCount, eThresh):
		retList = list()
		speciesSortedMap = self.sortedBySpecies(genesOnly= True)
		for s in speciesSortedMap:
			blastHit = BLASTPSearch.searchListBLASTpBySpecies(s, speciesSortedMap[s], returnCount, eThresh, printOut = False)
			if blastHit:
				retList.extend(blastHit)
			if len(blastHit) < 2 or blastHit[0][2] == "Not in database":
				print "Some motifs were not available in the database."
				#print blastHit[0]
				continue
		return retList

