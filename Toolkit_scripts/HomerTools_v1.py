#!/usr/bin/python
###########################################
#This is a simple tool to perform useful functions having to do with the HOMER motif recognition tool.
##########################################
from Bio import SeqIO
import re
import os
import argparse

START = ">"

#An object for reading in and manipulating homer data, but most especially for creating
#custom homer files
#Each HomerData object represents a single motif.
class HomerData:
	def __init__(self, seq):
		self.seqName = seq
		if seq[0] != START:
			self.seqName = START + seq
		self.name = "ID"
		self.thresh = 8
		self.matrix = list()
	
	def write(self, openFile):
		openFile.write(self.seqName + '\t' + self.name + '\t' + str(self.thresh) + '\n')
		self.writeMatrix(openFile)
	
	def writeMatrix(self, openFile):
		for row in self.matrix:
			for i in range(0, len(row)):
				if i == (len(row)-1):
					openFile.write(str(row[i]) + '\n')
				else:
					openFile.write(str(row[i]) + '\t')
					


###################Parse HTML from HOMER output################

#Helper function to write it out to a page.
def writeMotifMatchHTML(path):
	refMotifMap = dict()
	#writeOut= open(outputFile, 'w')
	for i in os.listdir(path):
		if i.endswith(".info.html"): 
			motifMatchHTMLExtraction(path + "/" + i, printOut = True)
	#print "Output written to", outputFile
	#writeOut.close()


#Main function: Extracts the list of matching motifs from the Homer output
#htmlFile- path to the file in
#@param outFile- the path to an output writing file if you want to do that.
#@return a list [seqeunce, list of each TFBS that matches that sequence]
def motifMatchHTMLExtraction(htmlFile, printOut = False):
	#The keys that we want
	refMotifDataList = list()
	htmlFileStream = open(htmlFile, 'r')
	currentLine = htmlFileStream.readline().strip()
	#get the second line
	currentLine = htmlFileStream.readline().strip() 
	#The Motif sequence is on the second line
	seq = extractMotifSeqHTML(currentLine)
	
	while currentLine != "":
		refMotifDataList.append(findMotifDataCoreHTML(htmlFileStream))
		currentLine = htmlFileStream.readline().strip()
	htmlFileStream.close()
	
	#Option to write these out.
	if printOut:
		#print refMotifDataList
		for subList in refMotifDataList:
			if subList is not None:
				print seq, '\t',
				for element in subList:
					if type(element) is str:
						print element, '\t',
				print ''
		return
		
	return [seq, filter(None, refMotifDataList)]  #filter out the leftover nonetype matches.
	

#Extracts the motif's sequence from the HTML file
#Useful as a key
def extractMotifSeqHTML(currLine):
	if "<H2>" in currLine:
		return re.findall("([A-Z\(\)]+)\s", currLine)[0]
		
#Aids motifMatchHTMLExtraction: actually parses for the information from the page that we care about
#@param fileStream:
#@param the motif sequence
#@return a list containing the [name, rank, match_score]
def findMotifDataCoreHTML(fileStream):
	regexOptions= ['<H4>\d+(M[\w\d\._-]+)<\/H4>', '<H4>\d+([\w\d\._-]+)<\/H4>','<H4>\d*([\w\d\._-]+)<\/H4>']
	ID_HEADER = "<H4>"
	ID = "MotifID"
	MATCH_RANK = "Match Rank:"
	SCORE = "Score:"
	OFFSET = "Offset:"
	motifData = list()
	currentLine = fileStream.readline().strip()
	while currentLine != "":
		if ID_HEADER in currentLine:
			header = re.search("<H4>\d+-([\w\d\._-]+)<\/H4>", currentLine)
			regexOption = 0
			#Make sure we get a valid response
			while header is None:
				header = re.search(regexOptions[regexOption], currentLine)
				regexOption += 1
			motifData.append(header.group(1))
		elif MATCH_RANK in currentLine:
			motifData.append(extractHTMLValue(currentLine, MATCH_RANK))
		elif SCORE in currentLine:
			motifData.append(extractHTMLValue(currentLine, SCORE))
		elif OFFSET in currentLine:  #After we get all the data we want
			if len(motifData) != 3:
				print "ERROR- didn't get everything"
			return motifData
		else:
			pass
		currentLine = fileStream.readline().strip()
	
			
#A specific tool to extract HTML values from the homer output html
#Specific to this
#@param currentLine of the HTML file
#@param key to look for and excise
#@return the value with that assignment
def extractHTMLValue(currentLine, key, group1 = False):
	regex = "<TD>([\w\s\.:]+)<\/TD"
	regexResult = regexExtraction(currentLine,regex,  groupNum = 2)
	if regexResult is not None:
		return regexResult
	else:
		return ""



#Simple function to get the name of a motif using regular expressions
#String to search, 
#@param regular expression to search with
def regexExtraction(strIn, regex, groupNum = 1):
	regex = re.findall(regex, strIn)
	if regex is not None:
		return regex[groupNum-1]
	else:
		return None


#####################################Creating custom HOMER files############################################

#The core function: itertes through a list of HOMER objects and writes them out in the appropriate format
#@param filePath: the file path you want to write to
#@param the list of homer objects
def writeHomerDB(filePath, homerList):
	HOMERDB = open(filePath, 'w')
	for item in homerList:
		item.write(HOMERDB)
	HOMERDB.close()

#A custom tool made for converting the downloadable library from PlantPan into a .motifs HOMER compatible database
#@param inFilePath- the PlantPan databse file
#@param outFilePath- where you wish to write the new HOMER compatible library to
#@customThresh-  the logOddsDetection threshold required for HOMER files.  For a Db, this value can be anywhere between 
##5 and 10 with little effect on outcome
def convertPPToHomerLibrary(inFilePath, outFilePath, customThresh = 8):
	PWSDB = open(inFilePath, 'r')
	motifList = list()
	currentLine = PWSDB.readline().strip()
	lineCount = 1
	regex = "\[([\w\t\.\^-]+)\]"
	while  currentLine != "":
		if currentLine[0] == START:
			#New 
			currHomer = HomerData(currentLine)
			#Get the next four lines
			a_line = extractMatrixLineValues(PWSDB.readline().strip(), regex)
			c_line = extractMatrixLineValues(PWSDB.readline().strip(), regex)
			g_line = extractMatrixLineValues(PWSDB.readline().strip(), regex)
			t_line = extractMatrixLineValues(PWSDB.readline().strip(), regex)
			currHomer.matrix = transformMatrix([a_line, c_line, g_line, t_line])
			currHomer.thresh = customThresh
			currHomer.name = str(lineCount) + "-" + currentLine[1:] 
			motifList.append(currHomer)
		currentLine = PWSDB.readline().strip()
		lineCount +=1
	writeHomerDB(outFilePath, motifList)
	print "Converted PlantPan database to HOMER format!  Check results."

#Converts CIS-BP format PWMs into a HOMER-compatible .motifs file
#@param inPath- the path to a directory containing all of the PWMs you wish to include in the HOMER file
##Note that this is the standard CIS-BP download format available on line
#@param outfilePath- where you would like to write the final .motifs database file to
#@param customThresh- the logOddsDetection threshold required for HOMER files.  For a Db, this value can be anywhere between 
##5 and 10 with little effect on outcome
def convertCISPBtoHomerLibrary(inPath, outFilePath, customThresh=8):
	motifList = list()
	regex = "\d\t([\w\t\.\^-]+)"
	counter = 0
	for i in os.listdir(inPath):
		currFile = open(inPath + "/" + i, 'r')
		currentLine = currFile.readline().strip()
		currMotifData = list()
		idName = regexExtraction(i, '(M\d{4}_\d\.\d\d).txt')
		currHomer = HomerData(idName)
		while currentLine != "":
			if "Pos" not in currentLine:
				currMotifData.append(extractMatrixLineValues(currentLine,regex))
			currentLine = currFile.readline().strip()
		#print "Finished file:", idName
		if len(currMotifData) == 0:
			continue			
		currHomer.matrix = currMotifData
		currHomer.thresh= customThresh
		currHomer.name = str(counter) + idName
		currFile.close()
		motifList.append(currHomer)
		counter += 1
	writeHomerDB(outFilePath, motifList)
	print "Converted CISBP database to HOMER format!  Check results."
		
				
				

#Simple a helper tool for safely extracting lines from a file
#@return a list containing each line in the file
#@param fileStream to extract from
def safeLineExtraction(fileStream, lineCount = 4):
	lineList = list()
	for i in range(0, lineCount):
		currLine = fileStream.readline().strip()
		if currLine == "":
			print "Empty line found"
		else:
			lineList.append(currLine)
	if len(lineList) == lineCount:
		return lineList

###########################################Matrix Operations#####################
#Takes a matrix- list of lists- and transforms it into HOMER accepted format
def transformMatrix(matrix):
	if matrix is None:
		print "Matrix is of noneType"
		return list()
		
	retMatrix = list()
	for i in range(0, len(matrix[0])):
		currVector = list() 
		for j in range(0, len(matrix)):
			currVector.append(matrix[j][i])
		retMatrix.append(currVector)
	return retMatrix
	
#Uses a custom regular expression to extract the matrix data from a given line
#@param line to search from
#@searchRegex to search using, can't be easily customized :(
#@delim: specify the delimiter between matrix values, default is tab
#@return: 
def extractMatrixLineValues(line, searchRegex, delim = '\t'):
	dataExtract = re.search(searchRegex, str(line))
	if dataExtract:
		return dataExtract.group(1).strip().split(delim)
	else:
		print "No data located", line
		return list()
	

###########################################################################################

################################Main function
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='HOMER tool with multiple functions.  htmlParser parses HOMER html output from compareMotifs.pl.  libConversion takes a database and converts it to .motifs file')
	parser.add_argument("-hs", "--htmlSource", help="Specify the path to the homer output of html files.")
	parser.add_argument("-o", "--output", default = "TF_match.tsv", help="Specify the comparison output file, default is \"TF_match.tsv\"")
	parser.add_argument("--htmlParser", action = "store_true", help = "Specify this if you want to use the html parsing tool and print an output")
	parser.add_argument("-l", "--libConversion",choices=["PP", "CISBP"], help = "Specify the type of library you wish to convert")
	parser.add_argument("-lp", "--libPath", help = "Specify the path to the PP or CISBP library you wish to convert")
	parser.add_argument("-log", "--logOddsDetection", type=int, default = 8, help = "Specify a log odds detection threshold to useif you are doing a library conversion. Should be between 5.0 and 10.0")
	args = parser.parse_args()
	
	if(args.libConversion == "PP"):
		if not args.libPath:
			convertPPToHomerLibrary("/home/likewise-open/ICE/aomdahl/Datasets/PlantPan/Transcription_factor_weight_matrix.txt", args.output, args.logOddsDetection)
		else:
			convertPPToHomerLibrary(args.libPath, args.output, args.logOddsDetection)
	if (args.libConversion == "CISBP"):
		if not args.libPath:
			convertCISPBtoHomerLibrary("/home/likewise-open/ICE/aomdahl/Datasets/CIS-BP/PWMs/pwms", args.output, args.logOddsDetection)
		else:
			convertCISPBtoHomerLibrary(args.libPath, args.output, args.logOddsDetection)

	if(args.htmlParser):
		writeMotifMatchHTML(args.htmlSource)
