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


####################################################################################################################
#This script simply merges output files written for 1 and 2kb sections and places them in a convenient folder.
#The files to merge must have the same name.
def fixSize(lengthIn, table):
	if len(table) > lengthIn:
			print "Weird case..."
			return
	while len(table) != lengthIn:
		table.append(" ")

if __name__ == '__main__':
	#Argument parsing
	parser = argparse.ArgumentParser(description='combines output files based on their directory')
	parser.add_argument("-f1", "--folder1", default = "./1kb/FinalReports", help="Folder 1on ")
	parser.add_argument("-1L", "--f1Length", default= "1kb", help = "Length of folder 1 data")
	parser.add_argument("-2L", "--f2Length", default= "2kb", help = "Length of folder 2 data")
	parser.add_argument("-f2", "--folder2", default = "./2kb/FinalReports", help="Second folder to merge")
	args = parser.parse_args()

	if not os.path.exists("CombinedFinalResults"):
		os.makedirs("CombinedFinalResults")
	
	try:
		for fPath in os.listdir(args.folder1):
			#Makes sure its a file
			if fPath.endswith(".tsv"):
				#Find the comparable file in the other folder
				print "Found first copy..."
				print fPath
				print os.path.join(args.folder2, fPath)
				if os.path.isfile(os.path.join(args.folder2, fPath)):
					print "Successfully located both copies of ", fPath
					#Parse them
					fileOne = open(os.path.join(args.folder1, fPath), 'r') 
					fileTwo = open(os.path.join(args.folder2, fPath), 'r')
					fileOut = open(os.path.join("CombinedFinalResults", fPath), 'w') 
					dataList = list()
				#print the top line
		
					currentLine = fileOne.readline().strip()
					headLength = len(currentLine.split('\t'))
					fileOut.write(currentLine + '\t' + "PromoterLength" + '\n')
					
					while currentLine != "":
						if "value" in currentLine or "Value" in currentLine:
							print currentLine
						#elif "NIAT" not in currentLine:
						#	pass
						else:
							data = currentLine.split('\t')
							fixSize(headLength, data)
							data.append(args.f1Length)						
							dataList.append(data)
							
						currentLine = fileOne.readline().strip()
						
					fileOne.close()
		
					currentLine = fileTwo.readline().strip()
					counter = 0
					while currentLine != "":
						if counter != 0:
							#if "NIAT" in currentLine:
							data = currentLine.split('\t')
							fixSize(headLength, data)
							data.append(args.f2Length)						
							dataList.append(data)
						counter += 1
						currentLine = fileTwo.readline().strip()
						
					fileTwo.close()
			
					sortedData = sorted(dataList, key=lambda x:x[0])
					
					for data in sortedData:
						fileOut.write(('\t').join(data) + '\n')
					fileOut.close()
				else:
					print "Unable to find this file in the neighboring folder.."
			else:
				print fPath
			

		#store the othe rlines
		#sort them by first key 
		#Write them out with the same file name in the current directory
	except OSError:
		print "Please make sure the given file paths are correct."
		sys.exit()
