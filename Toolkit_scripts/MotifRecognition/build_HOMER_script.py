#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#Takes a list of promoter sequences and organism genome 
#to generate a script to perform motif analysis
#
#Copyright 2016 aomdahl <aomdahl@ECO-12>
#
#

import stat
import sys, commands
import argparse
import os
sys.path.insert(0,"/home/likewise-open/ICE/aomdahl/Toolkit")
import GeneUtils as GU
import random

GENERAL_DIRECTORY = "MotifIdentification"
BOOTSTRAP_DIRECTORY = "PromoterBootstrapSequences"
PMOTIF_DIRECTORY =  "PutativeMotifs"
BOOTSTRAP_MOTIF_DIRECTORY = "BootstrapMotifResults" #Removed MotifIdentification/, add back in if you get an issue
MOVE_UP = ".."
HOMER_SCRIPT = "MotifIdentificationScript"



#This function creates a file with random boostrap promoter sequence ("false positives") from different genes
#For comparison and validation of actual data
#@pararm seqDict-- the list of promoter sequences in a dictionary
#@param cycleRepeat- number of times to repeat
#@param promtorRepeat- number of promoters to generate on each run
def randomBootstrapBuilder(posSeqDict, bgSeqDict, cycleRepeat):
    os.chdir(BOOTSTRAP_DIRECTORY)
    print "Constructing negative bootstrap dataset..."
    positivePromoters= set(posSeqDict.keys())
    #print "list length", len(positivePromoters)
    #print "set length,", len(posSeqDict.keys())
    for x in range(0,cycleRepeat):
        fileName = "round"+str(x)+".fa"
        output = open(fileName, 'w')
        
        for y in range(0,len(positivePromoters)):
            seqID = random.choice(bgSeqDict.keys())
            if seqID not in positivePromoters:
                output.write(">" + seqID + '\n')
                output.write(bgSeqDict[seqID]+ '\n')
            else:
                #print seqID + "will be omitted from bootstrap sequences (repeat or positive)"
                continue
            #usedPromoters.append(seqID)
        output.close()
        os.chmod("./"+fileName, stat.S_IRWXG)
        os.chmod("./"+fileName, stat.S_IRWXO)
        os.chmod("./"+fileName, stat.S_IRWXU)
    os.chdir(MOVE_UP)

#This actually writes the HOMER script
#@param proSeq is the list of promoters(a dictionary)
#@param dbPath is the path to the genome scaffold library
#@param mLength is the length of motif to look for
#@param processors is the number of processors to use
#@outpath is the place to put the output file
#Syntax is as follows:
    #"findMotifs.pl promoters fasta outputFile -fastaBg backgroundNoise -len lengthofMotif .... -p Number ofProcessors"
def scriptHOMER(writeFile, proSeq, dbPath, mLength, processors, outPath):
    cmd = "findMotifs.pl " + proSeq + " fasta " + outPath + " -fastaBg " + dbPath + " -len " + str(mLength)\
            + " -bits -basic -nogo -nocheck -p " + str(processors) + " -nofacts"
    writeFile.write(cmd + '\n')





#This function writes out the script to make the homer call
#proSeq: the path to the list of promoter sequences
#dbPath: the path to the plant's genome dbPath
#mLength: the length of motif to look for
def writeHOMERScript(proSeq, dbPath, mLength):
	scriptFile = open(HOMER_SCRIPT + str(mLength) + ".sh", 'w')
	#First, write the script to call the normal homer run.
	scriptHOMER(scriptFile, proSeq, dbPath, mLength, processors,PMOTIF_DIRECTORY+"/experimentalResults")
	#print "Final call", dbPath
	#Second, write the script lines to call the bootstrap runs
	try:
		bootstrapFiles = os.listdir(BOOTSTRAP_DIRECTORY)
		for boot in bootstrapFiles:
			roundName = boot.replace(".fa","")
			scriptHOMER(scriptFile, BOOTSTRAP_DIRECTORY + "/" + boot, dbPath, mLength, processors, 
				PMOTIF_DIRECTORY+"/"+ BOOTSTRAP_MOTIF_DIRECTORY+"/"+roundName)
	except OSError:
		scriptFile.close()
		return
		
	
	scriptFile.close()

#Tool to safely construct a new directory for this
#@param path
def makeDirectory(path):
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)

#After instantiating all of the command line arguments,
#build the random bootstrap samples and write the homer script

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description= "This script creates a script to perform a modified HOMER2 search, by using a provided gene database and promoter sequence list to find overrepressented motifs against a negative random sample from the plant genome. \n \t This also does a bootstrap comparison to confirm validity.")
	parser.add_argument("proSeq", help ="Promoter sequences, in a FASTA file")
	parser.add_argument("database", help = "Path to a promoters database, used for random background noise.  FASTA file")
	parser.add_argument("-r", "--runs", type = int, default = 100, help = "The number of random boostrap runs to perform, default is 100")
	parser.add_argument("-p", "--processors", help= "Number of processors to use, default is 2")
	parser.add_argument("-l", "--motifLength", help= "Length of the motif to search for, default is 8")
	parser.add_argument("-nb", "--noBoots", action="store_true", default=False, help = "Select this option to create scripts with no bootstrap runs")
	args = parser.parse_args() 
	#Create the output file for the program.
    
	global processors
	global motifLength
	global dbPath

	if args.processors:
		processors = args.processors
	else:
		processors = 2
    
	if args.motifLength:
		mLength = args.motifLength
		GENERAL_DIRECTORY = GENERAL_DIRECTORY + str(mLength)
	else:
		mLength = 8
        
	motifLength = mLength
	makeDirectory(GENERAL_DIRECTORY)
	makeDirectory("./"+GENERAL_DIRECTORY + "/"+ PMOTIF_DIRECTORY)
            
	global CORE_DIRECTORY
	CORE_DIRECTORY = os.path.abspath(os.getcwd())
       
	dbPath = os.path.abspath(args.database)
	proSeqPath = os.path.abspath(args.proSeq)

	makeDirectory(GENERAL_DIRECTORY + "/"+BOOTSTRAP_DIRECTORY)
	os.chdir(GENERAL_DIRECTORY)
	
	if not args.noBoots:
		randomBootstrapBuilder(GU.genomeScaffoldMap(open(proSeqPath,'r')), GU.genomeScaffoldMap(open(dbPath, 'r')), args.runs)
		
	writeHOMERScript(proSeqPath, dbPath, mLength)
