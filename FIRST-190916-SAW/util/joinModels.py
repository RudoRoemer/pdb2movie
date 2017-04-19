#!/usr/bin/python
import sys
import re

if (len(sys.argv) != 3) :
	print "Usage: " + sys.argv[0] + " inputSetName outputFilename "
	sys.exit()

inputFileName = sys.argv[1]
outputFileName = sys.argv[2]

inputList = open(inputFileName, 'r')

modelNumber = 1

outputFile = open(outputFileName, 'w')

for filename in inputList :
  filename = filename.strip()
  print filename
  inputFile = open(filename, 'r')

  outputFile.write("MODEL % 9d\n" % modelNumber)

  for line in inputFile :
    outputFile.write(line)

  inputFile.close()

  outputFile.write("ENDMDL\n")
    
  modelNumber += 1

inputList.close()
