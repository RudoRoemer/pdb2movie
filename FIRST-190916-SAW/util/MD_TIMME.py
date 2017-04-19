#!/usr/bin/python

import os
import sys

from optparse import OptionParser
optionParser = OptionParser()
optionParser.add_option("-f", "--frames", type="int", dest="numberOfFrames", help="The total number of frames output for TIMME analysis", default=50)
optionParser.add_option("-s", "--start", type="int", dest="startFrame", help="The first frame to output for TIMME analysis", default=0)
optionParser.add_option("-e", "--end", type="int", dest="endFrame", help="The number of frames before the last frame to output for TIMME analysis", default=0)

(options, args) = optionParser.parse_args()

VMD = "/home/smenor/bin/vmd"
FIRST = "/home/smenor/FIRST/src/FIRST"

numberOfFrames = options.numberOfFrames
startFrame = options.startFrame
endFrame = options.endFrame

VMDSCRIPTFILENAME = "vmdScript.tcl"

VMDScriptFile = open(VMDSCRIPTFILENAME, 'w')

VMDScriptFile.write("""
proc exportPDBfiles  {{outputFilePrefix prefix } {desiredNumberOfFrames 20} {startFrame 0} {endFrame 0}} {
  set mol 0
  set numberOfFrames [molinfo $mol get numframes]
  set selection [atomselect top all]

  set frameStep [expr (($numberOfFrames-$endFrame) - ($startFrame))/$desiredNumberOfFrames ]

  set multiModelPDBcontent ""
  set model 0

  for {set frame $startFrame} {$frame < $numberOfFrames - $endFrame} {incr frame $frameStep} {
    set multiModelPDBcontent [format "%sMODEL% 9d\n" $multiModelPDBcontent $model]
    
    $selection frame $frame
    
    set currentModelFileName [format "%s.%05d.pdb" $outputFilePrefix $model ]

    $selection writepdb $currentModelFileName

    # wait until the pdb file is written 
    while {True} {
      after idle [break]
    }

    set fd [open $currentModelFileName r]
    set modelContent [read $fd]
    close $fd

    set multiModelPDBcontent [format "%s\\n%s\\n" $multiModelPDBcontent $modelContent]

    set multiModelPDBcontent [format "%sTER\\nENDMDL\\n" $multiModelPDBcontent]

    exec rm -f $currentModelFileName
    incr model
  }

  set fd [open [format "%s.%d.models.pdb" $outputFilePrefix $desiredNumberOfFrames] w]
  puts $fd $multiModelPDBcontent
  close $fd
}
""")

extensions = []

outputFilename = args[-1]

for argument in args[:-1] :
  extension = argument.split(".")[-1]
  extensions.append(extension)
  
  if (extension == "pdb") :
    VMDScriptFile.write("mol addfile "+argument+"\n")
  
  if (extension == "nam") :
    namFile = open(argument, 'r')
    
    for line in namFile :
      filename = line.strip()
      
      VMDScriptFile.write("mol addfile "+filename+"\n")
  
  if (extension == "psf") :
    VMDScriptFile.write("mol addfile "+argument+"\n")
  
  if (extension == "top" or extension == "prmtop" or extension == "parm7") :
    VMDScriptFile.write("mol addfile "+argument+"\n")
  
  if (extension == "dcd") :
    VMDScriptFile.write("animate read dcd "+argument+" waitfor all\n")
     
  if (extension == "crd" or extension == "mdcrd") :
    VMDScriptFile.write("animate read crd "+argument+" waitfor all\n")

VMDScriptFile.write("exportPDBfiles %s %d %d %d \n" % (outputFilename, numberOfFrames, startFrame, endFrame))

VMDScriptFile.write("quit\n")

VMDScriptFile.close()

os.system(VMD + " -dispdev text -e " + VMDSCRIPTFILENAME)

print "run TIMME with the following: ", 
print FIRST, "-non -TIMME %s.%d.models.pdb" % (outputFilename, numberOfFrames)

