#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#USAGE: see -h for help options :)
#  

import string, sys, commands,os
import subprocess as sp
sys.path.insert(0,"/home/likewise-open/ICE/aomdahl/Toolkit")
import GeneUtils as GU
import argparse
import numpy
import errorReporting as ER
from sets import Set
from pMotif import pMotif
import cPickle as pickle
from MotifList import MotifList

#Map keys
CS = "CandidateSequence"
DEATS = "Details"
SOT = "ScoresOverThreshold"
WQ = "WithinQuantile"
TRUE = "True"
FALSE = "False"
EXACTS = "ExactMatchCount"
PVAL = "p-Value"
POS_FREQ = "PositiveFrequency"
NEG_FREQ= "NegativeBackgroundFrequency"
HFREQ_DIFF = "HomerFrequencyDifference"
PMATRIX= "ProbabilityMatrix"
LODT = "LogOddsDetectionThreshold"

#CompariMotif file keys/access
HEADER_START = "MotifFile"
EXACT_OVERLAP = "Exact Match"
EMPTY = "NONE"

#These are the dynamically identified indicies of the compariMotif file
#Output is variable: allows for flexibility for future versions
NAME_INDEX=-1
SCORE_INDEX=-1
MOTIF_INDEX=-1
DEATS_INDEX=-1
BOOT_NAME_INDEX =-1
BOOT_MOTIF = -1
OVERLAP_TYPE = -1

#Cutoff Thresholds
PERFECT_MATCH = 1
EXACT_MATCH_THRESHOLD = float(0.15) #If we have an exact match, its score must be below this threshold.
SCORE_THRESHOLD = float(0.5) #Threshold for a normalized similarity score (on a scale from 0 to 1)
PERCENTILE_THRESHOLD= float(0.05)	#Threshold for % of permitted overlaps with bootstrapped scores (i.e. 5% of 600 total bootstrap
									#sequences means a motif overlaps with 30 different bootstrap motifs with a score above the 
									#specific score threshold.
									#30.6 TEST CASE: 5% seems to be good-- keep samples, without isolating too many.
									
 
def run_command(cmd,ommit):
	if ommit:
		try: process = sp.Popen(cmd,shell=True)
		except: pass
		process.communicate("Y\n")
	else:
		try: process = sp.Popen(cmd,shell=True)
		except OSError,e: sys.exit("Error: Execution cmd failed")
		process.communicate("Y\n")
		
	return 0 


#Simple way to get the number of unique bootstraps from the file
#@param path: the path to the file to count the lines
def getEntryCount(path):
	val = sum(1 for lin in open(path, 'r'))
	return val


#Checks if our overlaps are within the specified percentile threshold
#@param length: number of bootstrap motifs with which the sequence has a score > the threshold score
#@param bootstrapEntryCount: total number of bootstrapEntries
#@return true if within the percentile, false if we are without- discard motif
def percentileCheck(length, bootstrapEntryCount):

	if float(length)> (PERCENTILE_THRESHOLD * bootstrapEntryCount):
	
		return False
	else:
		return True
		
#Does the current motif overlap with too many of the background motifs beyond the quantile
#threshold?
#@return True if so.
#@param motifMap the map of motifs
#@param motifName the motif in question
def exceedsPercentile(motifMap, motifName):
	if motifName in motifMap:
		if motifMap[motifName][SOT][WQ] == FALSE: 
			return True
	else:
		return False

#Safely increment the number of exact matches associated with a given motif
#@param motifKey the motif 
#@param motifMap the map of motifs
def incrementExactMatches(motifKey, motifMap):
	if EXACTS in motifMap[motifKey][SOT]:
		motifMap[motifKey][SOT][EXACTS] += 1
	else:
		 motifMap[motifKey][SOT][EXACTS] = 1

#Stores the motifs from the compariMotif file output into a map
#@param cmpareMotif_file: input file
#@bootstrapEntryCount: the number of bootstrap entries to compare against (for % calculation)
#cmRow stands for compariMotif file row
#Return the map of Candidate motifs
def buildMotifMap(cmpareMotif_file, bootstrapEntryCount, candidateMotifs):
	if len(candidateMotifs) == 0:
		return candidateMotifs
	exactMatches={}
	#Parse the file
	for line in open(cmpareMotif_file,"r"):#
		line = line.strip()
		cmRow = line.split("\t")
		if cmRow[0] == HEADER_START: #Check if we are on the top row of the file
			#Assign the appropriate access indices
			NAME_INDEX= GU.findSpecIndex(cmRow, "Name1")
			SCORE_INDEX = GU.findSpecIndex(cmRow, "Score")
			MOTIF_INDEX = GU.findSpecIndex(cmRow, "Motif1")
			DEATS_INDEX = GU.findSpecIndex(cmRow, "Desc1")
			BOOT_NAME_INDEX = GU.findSpecIndex(cmRow,"Name2")
			BOOT_MOTIF = GU.findSpecIndex(cmRow, "Motif2")
			OVERLAP_TYPE = GU.findSpecIndex(cmRow, "Sim1")
		
		else: #We are in the data of the file- not the header line
			motifName = cmRow[NAME_INDEX]
			motifSeq= cmRow[MOTIF_INDEX]
			#If there is a motif in the compariMotif output file that wasn't an original candidate, we have an issue
			if candidateMotifs[motifName][CS] != motifSeq:
				ER.storeError("TopMotif", "Motif sequences different:" + motifName, "Misaligned data files")
		
			#If the NORMALIZED score is greater than the score threshold, record it
			MAX_SCORE = len(cmRow[MOTIF_INDEX])
			MIN_SCORE = float(0)
			normScore = GU.normalizeScore(MIN_SCORE, MAX_SCORE, cmRow[SCORE_INDEX])
			
			#Store information regarding the scores exceeding the set threshold
			if float(normScore) >= SCORE_THRESHOLD:
				if cmRow[OVERLAP_TYPE] == EXACT_OVERLAP:
					incrementExactMatches(motifName, candidateMotifs)
				#Bootstrap name is the key referencing a table containing: [0] The normalized score, [1] motifSequence of Bootstrap, [2] the type of overlap
				candidateMotifs[motifName][SOT][cmRow[BOOT_NAME_INDEX]]= [normScore, cmRow[BOOT_MOTIF], cmRow[OVERLAP_TYPE]]

			else: #Score was lower than the cutoff
				if cmRow[OVERLAP_TYPE] == EXACT_OVERLAP: #Do we have an exact overlap?
					if  normScore >= EXACT_MATCH_THRESHOLD: #Does it exceed the "exact match" threshold condition?
						print "Exact match exceeding exact match threshold ", normScore
						candidateMotifs[motifName][SOT][cmRow[BOOT_NAME_INDEX]]= [normScore, cmRow[BOOT_MOTIF], cmRow[OVERLAP_TYPE]]
						print candidateMotifs[motifName][SOT][cmRow[BOOT_NAME_INDEX]]
			#If we have too much overlap with bootstrap sequences: i.e. motif has > 5% incidence with bootstrap sequences over threshold score
			if  not percentileCheck(len(candidateMotifs[motifName][SOT]), bootstrapEntryCount):
				candidateMotifs[motifName][SOT][WQ]=FALSE
			
	return candidateMotifs

#Generates a report in the event that motifs are found
def noMotifsReport():
	try:
		report = open("../FinalReports/ProcessReport.txt", 'w')
	except:
		report = open("ProcessReport.txt", 'w')
	report.write("No motifs passed filter thresholds.")
	report.write("Program terminated before further analysis.")
	report.write("Please see MotifIdentification*/PutativeMotifs/experimentalResults/homerMotifs.motifs* for possible motif candidates")
	report.close()





#Records the actual motifs that are candidates (pass all the threshold tests) in an output file
#@param outPath: output file, opened
#@param motifMap: the map of motifs
#@param motifDetails: the details to print alongside the 
#@param showAll- if we want to have access to all of the motifs, regardless of wether or not they got hit.
def writeTopMotifs(outPath, homerOutPath, motifMap, bootstrapEntryCount, showAll = True):
	
		
	#the output for Human viewing
	outputH = open(outPath, 'w')
	outputH.write("Name\tSequence\t%FreqInGeneSet\t%FreqInBackground\tFreqDiffHOMER\tExactMatchCount\tp-Value\tBootstrapOverlap\n")
	
	#The output file for the motif locating tool
	#We append here in the case of examining multiple motif lengths-- one super list with everything.
	outputSys = open(homerOutPath, 'a')
	
	#We want candidates to be printed in order of p-value
	#So we extract p-values, sort them by that, and then get the keys
	pValList=list()
	if len(motifMap) == 0:
		print "No candidate motifs identified for this run."
		print ".Motif file for conservation analysis written to ", homerOutPath
		outputH.close()
		outputSys.close()
		
	for candidate in motifMap:
		addList = [candidate, float(motifMap[candidate][DEATS][PVAL])]
		pValList.append(addList)

	highToLow = sorted(pValList, key=lambda x:x[1])

	motifMatchCount = 1
	for candidateTuple in highToLow:
		candidate=candidateTuple[0]
		if not showAll:
			if exceedsPercentile(motifMap, candidate):  #Don't print it
				continue
		else:
			#This output is exclusive to human viewing, and so can show anything.
			outputH.write(candidate + '\t')
			outputH.write(motifMap[candidate][CS]+'\t')
			outputH.write(motifMap[candidate][DEATS][POS_FREQ] + '\t')
			outputH.write(motifMap[candidate][DEATS][NEG_FREQ] + '\t')
			outputH.write(str(motifMap[candidate][DEATS][HFREQ_DIFF]) + '\t')
			printExactMatches(outputH, motifMap, candidate, bootstrapEntryCount)
			outputH.write(motifMap[candidate][DEATS][PVAL] + '\t')
			outputH.write(str(len(motifMap[candidate][SOT])/float(bootstrapEntryCount)) + '\t')
			outputH.write('\n')
			outputH.write(motifMap[candidate][DEATS][PMATRIX] + '\n')
			
			#This file can only show the eligible candidates- used for further analysis
			if not exceedsPercentile(motifMap, candidate):
				outputSys.write(">" + motifMap[candidate][CS] + '\t')
				outputSys.write(str(motifMatchCount) + "-" + motifMap[candidate][CS] + '\t')
				outputSys.write(str(motifMap[candidate][LODT]) + '\n')
				outputSys.write(motifMap[candidate][DEATS][PMATRIX] + '\n')
				motifMatchCount += 1
			
			#Store the motifs for serialization
			finalMotifList.addMotif(pMotif(motifMap[candidate]))

	print "Detailed output written to ", outPath
	print ".Motif file for conservation analysis written to ", homerOutPath
	outputH.close()
	outputSys.close()


#Small helper function to help with printout of exact matches
#@param output: output file, opened
#@param motifMap: the map of motifs
#@param motifKey: current key
def printExactMatches(output, motifMap, motifKey, bootstrapEntryCount):
	output.write(str(motifMap[motifKey][SOT][EXACTS]) + "/" + str(bootstrapEntryCount) + '\t')

#A simple function to store our candidate motifs in a map
#@param the path to the candidate motif file
#@return a map of name to sequence
def storeCandidates(candiPath):
	positives = open(candiPath, 'r')
	retData = dict()
	for line in positives:
		if EMPTY in line:
			return retData
		table = line.strip().split('\t')
		motifName = table[0]
		retData[motifName] = {}
		retData[motifName][SOT]={}
		retData[motifName][SOT][WQ]=TRUE #The score for a blank entry is not over threshold (SOT) yet
		retData[motifName][SOT][EXACTS] = 0 
		retData[motifName][DEATS]= {}
		#Assign the actual seqence
		retData[motifName][CS]=table[1]
	positives.close()
	return retData
		
#this function refers to the original homer file to extract the data
#we are interested in for later analysis
#@param filtCan: dictionary of filtered Candidate motifs
#@param homerPath: the path to the homer file
#@return the map containing the detail information.
def homerMetricsExtraction(filtCan, homerPath):
	hFile = open(homerPath, 'r')
	currentLine = hFile.readline().strip()
	if len(filtCan) == 0:
		return filtCan
	while currentLine != "":
		
		#If we are at a metric line
		if currentLine[0]== '>':
			for name in filtCan:
				if filtCan[name][CS] in currentLine:
					#Extract information
					homerHeader = currentLine.split('\t')
					metricsTable = (homerHeader)[5].split(',')
					#Assign the P-value
					filtCan[name][DEATS][PVAL] = GU.extractPValue(metricsTable[2])
					#Assign the background negative frequency
					filtCan[name][DEATS][NEG_FREQ] = GU.extractPercentage(metricsTable[1])
					#Assign the Positive Percentage
					filtCan[name][DEATS][POS_FREQ] = GU.extractPercentage(metricsTable[0])
					#Assign the Difference
					filtCan[name][DEATS][HFREQ_DIFF] = float(GU.extractPercentage(metricsTable[0])) - float(GU.extractPercentage(metricsTable[1]))
					
					#Extract the log odds detection threshold, used for making the output file for the next step.
					filtCan[name][LODT] = float(homerHeader[2])

					#extract the matrix
					matrix = ""
					for x in range(0, len(filtCan[name][CS])):
						currentLine = hFile.readline()
						matrix += currentLine
					#Assign the Matrix
					filtCan[name][DEATS][PMATRIX] = matrix
					
					 
		currentLine = hFile.readline().strip()
	hFile.close()
	return filtCan
			

#This tool examines motif overlap between the HOMER bootstrap entry and positive dataset
#This then uses that to calculate which motifs are significant and which are not
#All motifs and their calculated scores are viewable in "allMotifs.txt", listed with their % overlap with the bootstrap motifs
#The motifs making the cutoff (overlap with <=5%) are printed out into "topMatches.motif" in the super directory
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='This script examines motif comparison data and eliminates motifs with high similarity to bootstrap motifs.')
	parser.add_argument("input", help="the compariMotif produced file.")
	parser.add_argument("candidate", help="the path to the candidate motifs file- presto format")
	parser.add_argument("bootstrap", help = "the bootstrap motif file- presto format")
	parser.add_argument("unmatched", help ="Path to the unmatched motifs from the compariMotif output")
	parser.add_argument("homer", help="Path to the HOMER output file containing the motifs")
	parser.add_argument("-o", "--output", help="the output path to store produced files.")
	parser.add_argument("-c", "--cutoff", help ="Cutoff similarity score for motifs and bootstrapped values.  Must be between 0 and 1")
	parser.add_argument("-m", "--module", default = "M1", help ="The name of the module we are working with.  Default is M1")

		
	args = parser.parse_args()
	global outpath
	global candipath
	global bootstrapEntries
	global unmatchedPath
	global bootstrapEntryCount
	global finalMotifList
	
	if args.cutoff:
		cutoff = float(args.cutoff)
		if cutoff > 1 or cutoff < 0:
			print "Invalid cutoff entry: Must be between 0 and 1"
			exit()
		SCORE_THRESHOLD = float(args.cutoff)
	
	if args.input and args.candidate and args.bootstrap and args.unmatched and args.homer:		
		inputPath = os.path.abspath(args.input)
		bootPath = os.path.abspath(args.bootstrap)
		candiPath = os.path.abspath(args.candidate)
		unmatchedPath = os.path.abspath(args.unmatched)
		homerPath = os.path.abspath(args.homer)
		
		if args.output:
			outPath = os.path.abspath(args.output)
			if not os.path.exists(outPath):
				os.mkdir(outpath)
		else: 
			outPath = os.getcwd()+"/allMotifs.txt"
			homerOutPath = os.path.abspath(os.pardir) +"/topMatches.motif"  #pardir refers to the parent directory.

		finalMotifList = MotifList()
		##Assign the motif's module
		finalMotifList.moduleName = args.module
		
		bootstrapEntryCount = getEntryCount(bootPath)
		originalCandidates = storeCandidates(candiPath)
		candidateMap = buildMotifMap(inputPath, bootstrapEntryCount, originalCandidates)
		filteredCandidates = homerMetricsExtraction(candidateMap, homerPath)
		writeTopMotifs(outPath, homerOutPath, filteredCandidates, bootstrapEntryCount)
		
		#Serialize our map of motifs for later use.
		GU.updateSerializedMotifData("../v1" + GU.SER_MOTIF_MAP, finalMotifList)
		
	else:
		print "Please specify compariMotif input file, candidate motif input file (.presto) and bootstrap motif file (.presto)"
