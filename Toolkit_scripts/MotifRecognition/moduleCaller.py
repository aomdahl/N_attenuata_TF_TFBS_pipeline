




import subprocess
import re
import argparse
import sys

parser = argparse.ArgumentParser(description='Generates commands to call modules.  Simple helper function')
parser.add_argument("-tf", "--tfList", default = "/home/likewise-open/ICE/aomdahl/Datasets/TranscriptionFactors", help="Specify the path to the file containing the list of transcription factor genes and their families")
parser.add_argument("-ml", "--moduleList",help = "Select this option to search for just a specific gene")
args = parser.parse_args()

if not args.tfList or not args.moduleList:
	print "Please specify both the list of transcription factors and the list of modules"
	sys.exit()

moduleF = open(args.moduleList, 'r')
tfList = args.tfList
for line in moduleF:
	if "NIAT" in line:
		currentLine = line.strip()
		geneID = currentLine[:13]
		try:
			hit = subprocess.check_output(["grep", geneID, tfList])
		except subprocess.CalledProcessError:
			hit = ["","NA"]
			
		if "grey60" in currentLine:
				print "Found grey, ungrouped genes", geneID
				continue
		try:
			colorS = re.search("_([a-z]+)_module", currentLine)
			color = colorS.group(1)
			
		except:
			print "no hits for", geneID
			continue
		familyHit = hit.split(',')[1]
		print "$scriptPath $modPath/" + geneID+ "_" + color + "_module.txt "+ familyHit + "-"+geneID + " " + color

writeOut.close()
			

