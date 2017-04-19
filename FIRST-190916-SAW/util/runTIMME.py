#!/usr/bin/python

import sys
import re
import os

PATH_TO_FIRST = "/home/flexdev/FIRST_preview"
PATH_TO_UTIL = PATH_TO_FIRST + "/util/"

KILFOIL2PDB = PATH_TO_UTIL + "kilfoil2pdb.py"
JOINMODELS = PATH_TO_UTIL + "joinModels.py"
FIRST = PATH_TO_FIRST + "/src/FIRST"

UNZIP = "/home/flexdev/bin/unzip"
UNCOMPRESS = "/home/flexdev/bin/uncompress" 
BUNZIP2 = "/home/flexdev/bin/bunzip2"
GUNZIP = "/home/flexdev/bin/gunzip"

# FIXME - add error handling
if (len(sys.argv) < 2) :
  print "usage: "
  print sys.argv[0] + " inputFilename"
  sys.exit(1)

inputFilename = sys.argv[1]

pathRe = re.compile("\/")
splitInputPath = pathRe.split(inputFilename)
inputPath = "/".join(splitInputPath[:len(splitInputPath)-1])

inputFilename = splitInputPath[len(splitInputPath) - 1]

TMP_PATH = inputPath + "/.tmp" 

os.chdir(inputPath)

dotRe = re.compile ("\.")
splitFilename = dotRe.split(inputFilename)

extension = splitFilename[len(splitFilename) - 1]
arguments = ["-non", "-TIMME", "-flexweb"]

if (len(sys.argv) > 2) :
  arguments.extend(sys.argv[2:])

usingTempDirectory = False

while (extension in ["zip", "Z", "tar", "gz", "tgz", "bz2", "pdb"]) :
  if (extension in ["pdb"]) :
    arguments.append(inputFilename)
    print FIRST+ " " + str(arguments)

    os.execv(FIRST, arguments)
    sys.exit()
  
  if (extension in ["zip", "tar"]) :
    os.mkdir(TMP_PATH)
    os.chdir(TMP_PATH)
    os.system("ln ../" + inputFilename + " ./")
    usingTempDirectory = True
  
  if (extension == "zip") :
    os.system(UNZIP + " " + inputFilename)
    splitFilename = splitFilename[:len(splitFilename) - 1]
    inputFilename = ".".join(splitFilename)
    extension = splitFilename[len(splitFilename) - 1]
    if (extension == "pdb") :
      continue
 
  if (extension == "Z") :
    os.system(UNCOMPRESS + " " + inputFilename) 
    splitFilename = splitFilename[:len(splitFilename) - 1]
    inputFilename = ".".join(splitFilename)
    extension = splitFilename[len(splitFilename) - 1]
    continue
  
  if (extension == "bz2" ) :
    os.system(BUNZIP2 + " " + inputFilename)
    splitFilename = splitFilename[:len(splitFilename) - 1]
    inputFilename = ".".join(splitFilename)
    extension = splitFilename[len(splitFilename) - 1]
    continue

  if (extension == "gz") :
    os.system(GUNZIP + " " + inputFilename)
    splitFilename = splitFilename[:len(splitFilename) - 1]
    inputFilename = ".".join(splitFilename)
    extension = splitFilename[len(splitFilename) - 1] 
    continue
  
  if (extension == "tar") :
    os.system("tar -xf " + inputFilename)
    
  if (extension == "tgz") :
    os.system("tar -xzf " + inputFilename + "")
  
  prefix = ".".join(splitFilename[:len(splitFilename) - 1])
  namFileName = prefix+".nam"
  os.system("touch making.nam")
  os.system("find . | grep -E '.txt|.pdb'> "+namFileName)

  commonExtension = ""

  namFileList = open(namFileName)
  for fileName in namFileList :
    fileName = fileName.strip()
    splitFileName = dotRe.split(fileName)
    extension = splitFileName[len(splitFileName)-1]
    if (commonExtension == "") :
      commonExtension = extension
    else :
      if (commonExtension != extension) :
        print "extensions don't match"
        os.system("touch extensionsDontMatch")
        sys.exit()

  workingFilename = ""

  if (commonExtension == "pdb") :
    workingFilename = prefix + ".pdb"
    os.system(JOINMODELS + " " + namFileName + " " + workingFilename)
    arguments.append("-perResidue")
    arguments.append("-minimumRigidClusterSize5")

  if (commonExtension == "txt") :
    workingFilename = prefix + ".pdb"
    os.system(KILFOIL2PDB + " " + prefix+".nam " + workingFilename)
    arguments.append("-perRigidLabel")
    arguments.append("-minimumRigidClusterSize5")
    arguments.append("-stripeThickness3")
    arguments.append("-separationBetweenStripes0")
  
  os.system("cp " + workingFilename + " ../")
  if (usingTempDirectory) :
    os.chdir("../")

  arguments.append(workingFilename)

  print FIRST+ " " + str(arguments) 

  os.execv(FIRST, arguments)
  sys.exit()


