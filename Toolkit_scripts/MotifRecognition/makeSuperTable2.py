#!/usr/bin/env python
#Formatting for SLimSuite CompariMotif tool
#
#Copyright 2016 aomdahl <aomdahl@ECO-12>
#

import sys
import os
import re
import argparse
import subprocess
import glob

####################################################################################################################
#This script simply merges output files written for 1 and 2kb sections and places them in a convenient folder.
#The files to merge must have the same name.
def fixSize(lengthIn, table):
	if len(table) > lengthIn:
			print "Weird case..."
			return
	while len(table) != lengthIn:
		table.append(" ")

#Basically it could be any number of things.
#zf-HD-NIATv7_g37829
def getFamilyName(thing):
	#From NIAT to the end is the core gene
	core_gene = thing[-13:]
	family_name = thing[0:-14]
	return family_name, core_gene

	#AUX-IAA-NIATv7_g1692

def processLine(line, gene, family):
	t = line.split('\t')
	if len(t) == 11:
		return currentLine + "\t" + coreGene + '\t' + familyName
	else:
		counter = 1
		MLS = 5
		REGULATING_MODULE = 6
		if "." not in t[MLS] or "1" not in t[MLS]: #THIS MEANS THERE IS probably some kind of text in there.
			t.insert[MLS]
			
		


if __name__ == '__main__':
	#Argument parsing
	parser = argparse.ArgumentParser(description='combines output files based on their directory')
	parser.add_argument("folder1", help="Folder 1on ")
	args = parser.parse_args()

	currDir = args.folder1
	#print currDir
	
	if not os.path.exists(currDir):
		print "FAIL " + currDir
	
		sys.exit()

	for subdir, dirs, files in os.walk(currDir):
		for directory in dirs:
			#print directory
			#raw_input()
			fPath = directory
		#Makes sure its a file
		#print fPath
		#print os.listdir( currDir + "/" +  fPath)
			#print directory
			if len(directory) > 1 and "CombinedFinalResults" in os.listdir( currDir + "/" +  directory):
				for  dPath in os.listdir( currDir + "/" +  directory + "/CombinedFinalResults"):
					#print dPath
					#raw_input()
					if "pTF_library.tsv" not in os.listdir( currDir + "/" +  fPath + "/CombinedFinalResults"):
						print "ERROR IN:" + currDir + "/" +  directory + "/CombinedFinalResults"
						raw_input()
					
					if dPath.endswith("pTF_library.tsv"):
						ptfFile = open(currDir + fPath + "/CombinedFinalResults/pTF_library.tsv", 'r')	
			
						headerLine = ptfFile.readline().strip()
						currentLine = ptfFile.readline().strip()
						familyName, coreGene = getFamilyName(fPath)
						while currentLine != "":
							printLine = processLine(currentLine, coreGene, familyName)
							print printLine
							currentLine = ptfFile.readline().strip()
		#	 print len(dirs)
		break
			
