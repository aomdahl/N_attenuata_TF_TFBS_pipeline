#!/usr/bin/env python

#A tool to keep track of unexpected errors
import sys

errorMap = {}

#Store an error for update
#@param ToolName is a string
#@param Details also a string: tell us about the problem
#@param errorType: give us an identifier [RunTime|??]
def storeError(toolName, details, errorType):
	if toolName in errorMap:
		errorMap[toolName].append(errorType + ": " + details)
	else:
		errorMap[toolName] = [(errorType + ": " + details)]
	print "Error: ", toolName, details, errorType 

def printErrorReport():
    for key in errorMap:
        print key
        for entry in errorMap[key]:
            print errorMap[key]

def printErrorReport(outFile):
    for key in errorMap:
        outFile.write(key)
        for entry in errorMap[key]:
            outFile.write(errorMap[key])


#Checks to see if a file name matches what it is expected to.
def fileNameCheck(file_name, std_name, program_name):
	if std_name in file_name:
		return
	else:
		print program_name+ ": incorrect directory/file access \"" + file_name + "\""
		storeError(program_name, "FileName incorrect", "File Location Error")
		sys.exit()
