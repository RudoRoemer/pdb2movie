#!/usr/local/bin/python

from sets import Set
import sys
import re

# FIXME - add error handling
if (len(sys.argv) < 3) :
  print "usage: "
  print sys.argv[0] + " inputListFilename outputPDBfilename" 
  sys.exit(1)

inputListFilename = sys.argv[1]
outputPDBfilename = sys.argv[2]

fileList = open(inputListFilename, 'r')

commonSites = Set()

modelNumber = 1

for inputFileName in fileList :
  currentSites = Set()
  #outputFile.write("MODEL% 9d\n" %modelNumber)

  inputFileName = inputFileName.strip()

  inputFile = open(inputFileName, 'r')

  for line in inputFile :
    lineComponents = line.split()

    numericLineComponents = []
    for component in lineComponents :
      numericLineComponents.append(float(component))
    
    siteNumber = numericLineComponents[4]
    residueNumber = siteNumber

    currentSites.add(siteNumber)

    #outputFile.write("ATOM  % 5d  C   BLB % 5d    %8.3f%8.3f%8.3f\n" %(siteNumber, residueNumber, numericLineComponents[0], numericLineComponents[1], numericLineComponents[2]))
    
  #outputFile.write("TER\n")
  #outputFile.write("ENDMDL\n")

  modelNumber += 1
  if len(commonSites) == 0 :
    commonSites = currentSites
  
  else :
    commonSites = commonSites.intersection(currentSites)

fileList.close()

fileList = open(inputListFilename, 'r')

outputFile = open(outputPDBfilename, 'w')
modelNumber = 1
for inputFileName in fileList :
  outputFile.write("MODEL% 9d\n" %modelNumber)
  
  inputFileName = inputFileName.strip()
  inputFile = open(inputFileName, 'r')

  for line in inputFile: 
    lineComponents = line.split()
    numericLineComponents = []
    for component in lineComponents :
      numericLineComponents.append(float(component))
  
    siteNumber = numericLineComponents[4]
    if (siteNumber not in commonSites) :
      continue
    residueNumber = siteNumber

    outputFile.write("ATOM  % 5d  C   BLB % 5d    %8.3f%8.3f%8.3f\n" %(siteNumber, residueNumber, numericLineComponents[0], numericLineComponents[1], numericLineComponents[2]))

  outputFile.write("TER\n")
  outputFile.write("ENDMDL\n")
  
  modelNumber += 1

outputFile.close()
