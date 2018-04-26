
#
#Start: 3 June 2016
#aomdahl
#

import os, sys, commands
import argparse
from Bio import SeqIO
import re
import gffutils
sys.path.insert(0,"/home/likewise-open/ICE/aomdahl/Toolkit")
import errorReporting as ER
import GeneUtils as GT #stands for gene toolbox
#Data/constant fields we want accessible
TOOL = "PROMOTER_LOCATOR"
genomeDict= {}
UP_DEFAULT = 2000
DOWN_DEFAULT = 0
US = "UPSTREAM"
DS = "DOWNSTREAM"
SENSE = "+"
ANTISENSE= "-"

def getPromoter(geneAnnot, proDim):
    #calculate the range of the promoterI
    if genomeDict[geneName(geneAnnot)] is None:
        ER.storeError(TOOL, "Missing scaffold data", geneName(geneAnnot))
        return None
    if geneAnnot.strand == SENSE:
        #NORMAL PROCEDURE
        upper = geneAnnot.start - proDim[US]
        lower = geneAnnot.start + proDim[DS]
        return promoterExtraction(upper, lower, genomeDict[geneName(geneAnnot)], SENSE) 
    if geneAnnot.strand == ANTISENSE:
        #OTHER PROCEDURE
        upper = geneAnnot.end + proDim[US]
        lower = geneAnnot.end - proDim[DS]
        #LOWER AND UPPER AND SWITCHED HERE SINCE ITS THE ANTISENSE STRAND!
        promoterSequence = promoterExtraction(lower, upper, genomeDict[geneName(geneAnnot)], ANTISENSE) 
        #print promoterSequence
        reverseStrand = GT.flipMinusStrand(promoterSequence)
        #print reverseStrand
        return reverseStrand

def geneName(geneAnnot):

    #print geneAnnot.seqid
    return geneAnnot.seqid 
    

#This actually gets the promoter sequence of nucleotides
def promoterExtraction(upper, lower, sequence, strand):
    if(strand == SENSE):
        if upper < 0: upper = 0
        if lower > len(sequence):
            lower = len(sequence)-1
        return sequence[upper:lower].upper()
    elif strand == ANTISENSE:
        if upper > len(sequence):
            upper = len(sequence)-1
        if lower < 0:
            lower = 0
        return sequence[upper:lower].upper()
    else:    
        #print "Upper", upper
    #print "Lower", lower
    #print "Sequence for processint",sequence
    #print sequence[upper:lower]
        return sequence[upper:lower].upper()

def writePromoter(geneAnnot, outputFile, proDim):
    outputFile.write( ">"+ geneName(geneAnnot) + '\n')
    outputFile.write(getPromoter(geneAnnot, proDim)+ '\n')
    return

#Sets the promoter size based on the user specs
#Otherwise, it sets it to the default sizes

def setPromoterSize(argParser):
    seqDim = {}
    if argParser.upstream_distance: 
        seqDim[US] = int(argParser.upstream_distance)
    else:
        seqDim[US] = UP_DEFAULT
    if argParser.downstream_distance:
        seqDim[DS] = int(argParser.downstream_distance)
    else:
        seqDim[DS]= DOWN_DEFAULT
    return seqDim

#A simple helper tool to give an appropriate output file name
def outputFile(args):
    if args.outputFile:
        return open(args.outputFile, 'w')
    else:
        fileNameStart = re.match( "[\w\.]*$", args.sequence)
        if fileNameStart is None:
            return open("PROMOTERS.fa", 'w')
        else:
            return open((fileNameStart + "PROMOTERS"), 'w')
        return open(fileName, "w")



#The main function
    #Inputs: Genome sequence data (.fa)  and corresponding annotations(.gff)
    #OutputL .fa fasta file with scaffold identifiers and promoter sequences
    #Promoter sequences are set to capture 1000bp upstream of the transcription start site and 
    #200 bp downstream of transcription start site

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description="This script will extract the promoter seqeuences from a gene list provided by the user. \n Note that the promoter sequece is predicted based on the start codon, since there is innaccuracyin using 5' UTR as a prediction tool (??)")
    parser.add_argument("sequence", help="The path to the seqeunce FASTA file")
    parser.add_argument("annotations", help="The path to the annotations .gff file")
    parser.add_argument("-u","--upstream_distance", help="Number of base pairs upstream of start codon to include in promoter sequence")
    parser.add_argument("-d","--downstream_distance", help = "Number of base pairs downstream of start codon to include in promoter sequence")
    parser.add_argument("-o","--outputFile", help= "Path of destination file to write to. Should be FASTA format.")
    args = parser.parse_args()
   
    #First, we want to get in the genome file so we have access to the scaffold names and such
    genomeDict = GT.genomeScaffoldMap(open(args.sequence)) 
        #genomeDict: a complete map of scaffold names to genetic DNA sequences    
    
    #Here, we extract all of the relevant annotation data, using the gff utils database tool
    if args.annotations:
        global annotationsDB
        annotationsDB = gffutils.create_db(args.annotations,dbfn= "annotDB", force=True)
        promoterDimensions=setPromoterSize(args)
        
        output = outputFile(args)
        for utr in annotationsDB.features_of_type('CDS', order_by='start'):
              
            writePromoter(utr, output, promoterDimensions)
    
    output.close()
    
