#!/usr/bin/env python
###########################################
#This is a simple tool to perform basic tasks on microArray data
##########################################
from Bio import SeqIO
import re, sys
import argparse
import glob
#import errorReporting as ER
import subprocess
import os
import numpy as np
from decimal import * #This allows us to get the accuracy we need

HEADER = "#HEADER#"

#Loads the mapping of ProbeIDs to GeneIDs into memory for quick access
#@param mapPath of the gene mapping
#@return retMap
def loadProbeIDMapping(mapPath):
	mapping = open(mapPath, 'r')
	retMap = {}
	for line in mapping:
		table = line.strip().split('\t')
		if len(table) < 1:
			continue
		try:
			retMap[table[0].replace("\"", "")] = table[1].replace("\"", "")
		except IndexError:
			print "Error at", line
			
	mapping.close()
	return retMap
	
def loadMicroArrayData(dataPath):
	mapping = open(dataPath, 'r')
	retMap = {}
	lineCounter = 0
	for line in mapping:
		table = line.strip().split('\t')
		if len(table) < 1:
			continue
		if lineCounter == 0:
			retMap[HEADER] = line.strip()
			lineCounter += 1
			continue
		try:
			retMap[table[0].replace("\"", "")] = [Decimal(i) for i in table[1:]]  #Need at least 14 dps?
		except IndexError:
			print "Error at", line
		lineCounter += 1
	mapping.close()
	return retMap


#Gets the new map of results.
#@param probeIDMap-- the mapping of probeIDs to the pertinent gene
#@param microArray -- the microArray results in a dictionary
def mapOutFinalArray(probeIDMap, microArray):
	retMap = dict()
	print "probeIDMap size", len(probeIDMap)
	print "MicroArray size", len(microArray)
	for probeID in microArray:

		if probeID in probeIDMap:
			if probeIDMap[probeID] in retMap:  #Special case where there are duplicates- a gene has already been assigne dits data
				#print "Comparing two options..."
				#raw_input(probeIDMap[probeID])
				#Calculate the average of the data for each...
				avgCur = sum(retMap[probeIDMap[probeID]])/(len(retMap[probeIDMap[probeID]]))
				avgNew= sum(microArray[probeID])/(len(microArray[probeID]))
				#If the average expression of the new one is bigger, keep it
				if avgNew > avgCur:
					retMap[probeIDMap[probeID]] = microArray[probeID]
				continue
				#Otherwise, don't change the assignment
			else:
				#print probeID
				#print "Adding in..."
				retMap[probeIDMap[probeID]] = microArray[probeID]
		else:
			print "ProbeID" , probeID, "cannot be found in the probeIDMap"
			
			continue
	print probeIDMap.keys()
	retMap[HEADER] = microArray[HEADER]
	return retMap

def writeOutFinalArray(writeMap, outfile):
	outStream = open(outfile, 'w')
	outStream.write(writeMap[HEADER] + '\n')
	for geneId in writeMap:
		outStream.write(geneId.replace("\"", "") + '\t')
		outStream.write(('\t').join([str(i) for i in writeMap[geneId]]) + '\n')
	outStream.close()
	
################################Main
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='HOMER tool with multiple functions.  htmlParser parses HOMER html output from compareMotifs.pl.  libConversion takes a database and converts it to .motifs file')
	parser.add_argument("-m", "--mapping", help="Path of mapping of probe IDs to geneIDs ")
	parser.add_argument("-d", "--data",  help="Actual path to micro array data")
	parser.add_argument("-o", "--output", help="Output file")
	args = parser.parse_args()
	if (args.mapping is None or args.data is None or args.output is None):
		print "Please specify all the necessary files.  See -h for help"
		sys.exit()
	probeIDs = loadProbeIDMapping(args.mapping)
	microArray = loadMicroArrayData(args.data)
	writeOutFinalArray(mapOutFinalArray(probeIDs, microArray), args.output)
	
