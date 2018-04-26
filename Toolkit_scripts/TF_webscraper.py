#!/usr/bin/env python
###########################################
#This is a simple tool to perform basic, rendundant functions used throughout the motif prediction toolkit
##########################################
import re
import sys
import os
import re
import argparse
import GeneUtils as GU
import urllib2
import xml.etree.ElementTree as ET
import time

global scrapeData
scrapeData = list()
##Key headers for output written to database:
ID= "TF_ID"
MOTIF= "Motif_ID"
TF = "TF_Name"
SPECIES = "TF_Species"
FAMILY = "Family_Name"
DEATS="Details"
nastyRegex= ">([\w\s\.!\\\/;_:#&%@=+()\[\]{}\*\|\',\^\"~,-]+[^<>\\t])<"

def extractMotifID(strIn):
	regex = re.search('(TF_motif_seq_\d*)', strIn)
	if regex is not None:
		return regex.group(1)
	else:
		return None

def extractMotifID(strIn, whole = False):
	regex = re.search('(\w+).xml', strIn)
	if regex is not None:
		if not whole:
			return regex.group(1)
		else:
			return regex.group(0)
	else:
		return None 

			
def webScrapePP(motifIdDb, output):
	##JUSt start by downloading all the web files.
	motifIDFile = open(motifIdDb, 'r')
	for line in motifIDFile:
		line = line.strip()
		motifID=""
		if line[0] == ">":
			motifID = extractMotifID(line)
			if motifID is not None:
				#do the onlnie quiery	
				print "Query:", "http://plantpan2.itps.ncku.edu.tw/TFBSinfo.php?matrix=" + motifID
				query = urllib2.urlopen("http://plantpan2.itps.ncku.edu.tw/TFBSinfo.php?matrix=" + motifID)
				data= query.read()
				outFile = open(output+"/" + motifID +  ".xml", 'w')
				outFile.write(data)
				outFile.close()
				time.sleep(1)
			if motifID is None:
				print "Un matching header:"
				print line
	
	##Cool, they are all downloaded.  Now what?
def parseXMLMatrixID(fileIn):
	xmlFile = open(fileIn, 'r')
	for line in xmlFile:
		line = line.strip()
		#horizontalTableParser(line, extractMotifID(fileIn)):
		if "GeneID=" in line:
			tableLines = line.split("</tr><tr>")
			#print tableLines
			for line in tableLines:
				#Iterate through matches in the line:
				resultIter = re.finditer(nastyRegex,line)
				captureList = list()
				while True:
					try:
						captureList.append(resultIter.next().group(1))
					except StopIteration:
						break
				newEntry = dict()
				if len(captureList) > 2:
					newEntry[ID] = captureList[0]
					newEntry[MOTIF]=extractMotifID(fileIn)
					newEntry[DEATS] = ""
					if len(captureList) == 3:
						newEntry[TF] = ""
						newEntry[SPECIES]= captureList[2]
						newEntry[FAMILY]=captureList[1]
					else:
						newEntry[TF] = captureList[1]
						newEntry[SPECIES]= captureList[3]
						newEntry[FAMILY]=captureList[2]
					scrapeData.append(newEntry)
				elif len(captureList) == 2 or len(captureList) == 1:
					print "Unusual case..."
					print extractMotifID(fileIn)
					raw_input("What's going on?")
					print captureList

	#print scrapeData
def parseXMLMotifSeqID(fileIn):
	fileLength = getLineCount(fileIn)
	xmlFile = open(fileIn, 'r')
	currentLine = xmlFile.readline().strip()
	newEntry = dict()
	matrixID = extractMotifID(fileIn)
	for i in range(0,fileLength-1): 
		motifSeqIDDataUpdate("Matrix ID", currentLine, xmlFile, newEntry, MOTIF)
		motifSeqIDDataUpdate("Matrix Name", currentLine, xmlFile, newEntry, TF)
		motifSeqIDDataUpdate("Species", currentLine, xmlFile, newEntry, SPECIES)
		motifSeqIDDataUpdate("Description", currentLine, xmlFile, newEntry, DEATS)
		newEntry[ID] = "NA"
		newEntry[FAMILY] = ""
		#Special case- for 136 of them...
		if "Correspounding TF" in currentLine:
			currentLine = xmlFile.readline().strip()
			currentLine = xmlFile.readline().strip()
			currentLine = xmlFile.readline().strip()
			currentLine = xmlFile.readline().strip()
			currentLine = xmlFile.readline().strip()
			horizontalTableParser(currentLine, matrixID)
		
		currentLine = xmlFile.readline().strip()
	
	#print "Appending the following entry..."
	#print newEntry
	scrapeData.append(newEntry)

	xmlFile.close()
	
	
	
#Current line we are on
#Matrix ID-- the id of the motif- file name.
def horizontalTableParser(line, matrixID):
	if "GeneID=" in line:
		tableLines = line.split("</tr><tr>")
		#print tableLines
		for line in tableLines:
			#Iterate through matches in the line:
			resultIter = re.finditer(nastyRegex,line)
			captureList = list()
			while True:
				try:
					captureList.append(resultIter.next().group(1))
				except StopIteration:
					break
			newEntry = dict()
			if len(captureList) > 2:
				newEntry[ID] = captureList[0]
				newEntry[MOTIF]=matrixID
				newEntry[DEATS] = ""
				if len(captureList) == 3:
					newEntry[TF] = ""
					newEntry[SPECIES]= captureList[2]
					newEntry[FAMILY]=captureList[1]
				else:
					newEntry[TF] = captureList[1]
					newEntry[SPECIES]= captureList[3]
					newEntry[FAMILY]=captureList[2]
				scrapeData.append(newEntry)
			elif len(captureList) == 2 or len(captureList) == 1:
				print "Unusual case..."
				print extractMotifID(fileIn)
				raw_input("What's going on?")
				print captureList

	
		
#Returns the assignment of what you want
def motifSeqIDDataUpdate(keyword, currLine, xmlFile, addMap, mapKey):
	if mapKey not in addMap:
		addMap[mapKey] = ""
	if keyword in currLine:
		dataLine = xmlFile.readline().strip()
		regexResult = re.search(nastyRegex, dataLine)
		if regexResult:
 			addMap[mapKey]= regexResult.group(1)
	return		

#Get the number of lines in teh file... not sure why we get the weird results.
def getLineCount(fileIn):
	xmlFile = open(fileIn, 'r')
	counter = 0
	for line in xmlFile:
		counter += 1
	xmlFile.close()
	return counter


#Takesour map of online scraped data and then viola 
#Writes it out to a file matching cis-bp output.
def writePPScrapeToFile(outputSpot):
	outputFile = open(outputSpot, 'w')
	#First print the headers
	outputFile.write("TF_ID\tMotif_ID\tTF_Name\tTF_Species\tFamily_Name\tDetails\n")
	for ex in scrapeData:
		outputFile.write(ex[ID] + '\t')
		outputFile.write(ex[MOTIF] + '\t')
		outputFile.write(ex[TF] + '\t')
		outputFile.write(ex[SPECIES] + '\t')
		outputFile.write(ex[FAMILY] + '\t')
		outputFile.write(ex[DEATS] + '\n')
	outputFile.close()






if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='This is a custom webscraping tool.')
	parser.add_argument("-i", "--motifIdDb", help="Path to the list of motif ids")
	parser.add_argument("-w", "--scrapedData", help="XML webscraper Output path")
	args = parser.parse_args()
	
	#webScrapePP(args.motifIdDB, args.output)
	if(args.scrapedData):
		for i in os.listdir(os.path.abspath(args.scrapedData)):
			idName = extractMotifID(i, whole=True)
			print "Current file:", idName
			if "TF_motif_seq_" in idName:
				parseXMLMotifSeqID(args.scrapedData + "/"+idName)
			elif "TFmatrix" in idName:
				parseXMLMatrixID(args.scrapedData + "/"+idName)
			else:
				print "Unknown file encountered:", idName
			 
	else:
		parseXMLMatrixID("/home/likewise-open/ICE/aomdahl/Datasets/PlantPan/WebData/TFmatrixID_0328.xml")
	#Parse the other kind of file
	
	#Write em up.
	writePPScrapeToFile("/home/likewise-open/ICE/aomdahl/Datasets/PlantPan/PlantPanMotifDB.tsv")

