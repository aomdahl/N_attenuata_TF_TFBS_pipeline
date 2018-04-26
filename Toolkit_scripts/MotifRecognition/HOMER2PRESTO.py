#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#Extracts motifs from HOMER output and puts into PRESTO format
#Formatting for SLimSuite CompariMotif tool
#
#Copyright 2016 aomdahl <aomdahl@ECO-12>
#
#
import os
import math
import re
import numpy
import argparse
import sys
sys.path.insert(0,"/home/likewise-open/ICE/aomdahl/Toolkit")

import GeneUtils as GU
import errorReporting as ER


BOOTSTRAP_FREQ_CUTOFF = 0.0 #The cutoff for motif frequency (matches) when running the bootstrap sequences-A.

MOTIF_LENGTH = 8 #Default value, is assigned as the files get parsed.

HOMER_COLUMN_COUNT = 7
HOMER_FILE_NAME= "homerMotifs"
EMPTY = "NONE"

#Checks the p-value of the sequence to make sure it meets the 
#specified threshold
#@param dataMetrics: string of relevant data from HOMER output
#@return -1 if it doesn't, pvalueF if it does
def checkPValue(dataMetrics, cutoff):
	table = dataMetrics.split(',')
	pvalueS = (table[2])[2:]

	pvalueF = float(pvalueS)
	if pvalueF > cutoff:
		return -1
	else:
		return pvalueF

#Makes sure the sequence appears in at least 20% of sequences
#9.6.2016: Removed check for background 50% threshold (i.e. don't include samples that are in more than
# 50% of backround data). This is the purpose of the bootstrap test, so the 50% cutoff is both arbitrary and redundant
#@param value: the string of motif detail information from the HOMER file
def setFilter(value):
	table = value.split(',')
	inSequences = table[0]
	#Use regular expressions to extract the percentages
	matchesPercent = percentageMatch(value)
	#Check to make sure the motif isn't just noise and appears with justifiable frequency.
	if matchesPercent is not None:
		if bootstrapMode:
			if float(matchesPercent) < BOOTSTRAP_FREQ_CUTOFF:
				return False
			return True
		else:
			if float(matchesPercent) < args.posFreq:
				return False
				
	else:
		ER.storeError("Homer2Presto", "Unable to extract data", "FileType")
		sys.exit()
	return True
		
#Helper method to quickly get the percentage of matches
#@param the percentage data string
#@return the percentage
def percentageMatch(dataMetrics):
	table = dataMetrics.split(',')
	inSequences = table[0]
	matchesPercent = re.search('[^\s]\((\d*\.\d*)%\)',inSequences)
	return matchesPercent.group(1)


##This loops through the HOMER output file
#@param fileStream of the HOMER output file
#@param sequenceList the map of motifs so far
def extractMotifs(fileStream, sequenceList):
	currentLine = fileStream.readline().strip()
	if len(currentLine.split('\t')) < HOMER_COLUMN_COUNT:
		ER.storeError("Homer2Presto", "Wrong number of columns in input file", "FileType")
		sys.exit()
	while currentLine != "":
		if currentLine[0] == '>':
			headers = currentLine.split('\t')
			
			#Extract the motif from the HOMER file
			motifSearch = (headers[1])
			motifRegex = (re.search('-[^A-Z]*(\w*)', motifSearch))
			motif = motifRegex.group(1)
			MOTIF_LENGTH= len(motif)
			
			#Extract the pValue from the HOMER file
			dataMetrics = headers[5]
			pValue = checkPValue(dataMetrics, args.cutoff)
		
			if pValue != -1 and setFilter(dataMetrics):
				safeAdd(sequenceList, motif, dataMetrics, pValue)

			
		currentLine = fileStream.readline().strip()


#A simple helper function to add quickly and securely to our map of sequences
#@param seqMap to add to
#@param motif to add
#@param dataMetrics: the percentage of sequences the motifs appears in
#Data added has the following structure [motifSequence, percentFreqency, Pvalue]
def safeAdd(seqList, motif, dataMetrics, pValue):
	percent = percentageMatch(dataMetrics)
	seqList.append([motif,percent, pValue])
	#July27 mods
	#if motif not in seqMap:
	#	seqMap[motif] = []
	#(seqMap[motif]).append(percent)


#Writes out the list of identified motifs in PRESTO format
#@param fileStream to write to
#@param seqMap the Map of motifs and their appearance %
#@param annotation- The name labels to accompany the motifs. Numbered from 1
def writeMapPRESTO(fileStream, seqList, annotation):
	if annotation:
		name = annotation
	else:
		name = "HomerMotifsNIATT"
	nameCounter = 1
	if len(seqList) == 0:
		fileStream.write(EMPTY + '\t' + EMPTY + '\t' + EMPTY + '\t')
		fileStream.close()
		return
		
	if not seqList:
		print "Empty list sent to presto creation"
		fileStream.close()
		return
	for motifArray in seqList:
		#[motifSequence, percentFreqency, Pvalue] is structure of motifArray
		fileStream.write(name + str(nameCounter)+ '\t')
		fileStream.write(motifArray[0] + '\t')
		#floatTable = [float(i) for i in seqMap[key]]
		#percent = numpy.mean(floatTable)
		fileStream.write("#" + str(motifArray[1]) + "% p-val:" +str(motifArray[2]) +  "\n")
		nameCounter += 1
	fileStream.close()




#This takes either a single HOMER output file or a tree of folders containg HOMER files and place them into one single array of motifs
#A motif is only kept if it appears in at least 30% of the promoter sequences
#July 27-- Modified this to use a list instead of a map, because we want all the unique ones
if __name__ == '__main__':
	#Argument parsing
	parser = argparse.ArgumentParser(description='Thi script extracts the relevant motif information from a HOMER output file (homerMotifs.all.motifs) \n'+
													'and is designed for use with a single positive sample and false positive boostrap samples')
	parser.add_argument("-b", "--bootstrap", help="The path to the folder containing the boostrap false positive samples")
	parser.add_argument("-i", "--inputFile", help="A single file to extract the relevant motifs from-- the ''positive'' dataset")
	parser.add_argument("-o", "--output", help="The output file to store list of candidate motifs.")
	parser.add_argument("-c", "--cutoff", type = float, default = 1e-5, help="The statistical significance threshold.  Default is 1e-5.")
	parser.add_argument("--posFreq", type = float, default = 20.0, help="Motif frequency threshold in positive dataset. Default is 20 percent.")
	parser.add_argument("-a", "--annotation", help= "Informative tag name to label the promoter sequnences [Bootstrap|CoExpressionGroupName]")
	args = parser.parse_args()
	if args.bootstrap and args.inputFile:
		print "Please specify either a directory or input file, not both"
		exit()

		
	#Keep track of which mode we are in-- bootstrap analysis or true positive promoter motif analysis
	global bootstrapMode
	sequenceList = list()
	#If we are examining the bootstrap data (a tree)
	if (args.bootstrap):
		
		print "Condensing bootstrap generated sequences into PRESTO format..."
		bootstrapMode = True
		bootPath = os.path.abspath(args.bootstrap)
		for root, dirs, files in os.walk(bootPath):
			for fileName in files:
				if fileName is not None and HOMER_FILE_NAME in fileName:
					
						currentFile = open(os.path.join(root,fileName), 'r')
						extractMotifs(currentFile, sequenceList)
						currentFile.close()
				
				
	#We are extracting from only one file-- the true positive data
	elif (args.inputFile):
		print "Condensing putative motifs into PRESTO format..."
		ER.fileNameCheck(args.inputFile, HOMER_FILE_NAME, "HOMER2PRESTO")
		bootstrapMode = False
		positivePath = os.path.abspath(args.inputFile)
		currentFile = open(positivePath, 'r')
		extractMotifs(currentFile, sequenceList)

	else:
		print "Missing -b|-f.  Please try again."
		exit()
		
	#Write the finalized motifs to the output file in PRESTO format.
	if (args.output):
		writeMapPRESTO(open(args.output, "w"), sequenceList, args.annotation)
	else:
		writeMapPRESTO(open("motifs.out.presto", "w"), sequenceList, args.annotation)
		

	
