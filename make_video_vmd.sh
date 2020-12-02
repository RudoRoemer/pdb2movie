#!/bin/bash

#check that five arguments have been given
if [ "$#" -ne 5 ]
then
	echo "usage: $0 video_width video_height frame_rate path_with_pdb_files destination_file_name_and_path"
	exit -1
fi

#arguments:
#PDB_PATH is where to find the PDBs which must begin 'tmp_froda_'.
#DEST_FILE is the name (and location) of the video to generate
VIDEO_WIDTH=$1
VIDEO_HEIGHT=$2
FPS=$3
PDB_PATH=$4
DEST_FILE=$5

TMP_PATH=$PDB_PATH/vmd_temp

#calculate how many seconds each frame lasts
FRAME_LENGTH=$(awk -v fps=$FPS 'BEGIN { print 1.0 / fps }')

#make a temporary folder for the temporary files
rm -rf $TMP_PATH
mkdir $TMP_PATH

#file to give vmd to generate all necessary frames (in tga format)
CMD_FILE=$TMP_PATH/vmd.cmd

PDB_FILES=$PDB_PATH/tmp_froda_*.pdb

#function allowing vmd zooming
echo 'proc AutoScaleAllVisible {{zoom_factor 1}} {
  set prog 1
  set me [lindex [info level [info level]] 0]
  #----
  set minX "false"
  set minY "false"
  set maxX "false"
  set maxY "false"
  set zoom [expr 1.8 * $zoom_factor]
 
  #compute bb for all visible stuff
  set mols [molinfo list]
  foreach molid $mols {
    if {[molinfo $molid get displayed] && [molinfo $molid get numframes] > 0} {
      set num_reps [molinfo $molid get numreps]
      if {$num_reps > 0} {
     set seltext ""
     for {set i 0} {$i<$num_reps} {incr i} {
       if {[mol showrep $molid $i]} {
         if {[string length $seltext] > 0} { set seltext "$seltext or " }
         set temp [molinfo $molid get "{selection $i}"]
         set seltext "${seltext}(${temp})"
       }
     }
      }
      if {[string length $seltext] > 0} {
     set sel [atomselect $molid ($seltext)]
     set mm [measure minmax $sel]
     $sel delete
     set minXtemp [lindex [lindex $mm 0] 0]
     set minYtemp [lindex [lindex $mm 0] 1]
     set maxXtemp [lindex [lindex $mm 1] 0]
     set maxYtemp [lindex [lindex $mm 1] 1]
     set minX [expr $minXtemp < $minX || $minX == "false" ? $minXtemp : $minX]
     set minY [expr $minYtemp < $minY || $minY == "false" ? $minYtemp : $minY]
     set maxX [expr $maxXtemp > $maxX || $maxX == "false" ? $maxXtemp : $maxX]
     set maxY [expr $maxYtemp > $maxY || $maxY == "false" ? $maxYtemp : $maxY]
      }
    }
  }
  if {$minX != "false"} {#true for 1 true for all
    set rangeX [expr $maxX - $minX]
    set rangeY [expr $maxY - $minY]
    set maxrange [expr max($rangeX,$rangeY)]
    set target [expr $zoom/$maxrange]
    set cscale [lindex [lindex [lindex [molinfo top get scale_matrix] 0] 0] 0]
    set nscale [expr $target + (($cscale - $target) * pow(abs($prog-1.),$prog))]
    eval "scale to $nscale"
  } else {
      puts "$me: nothing seems to be visible!"
  }
}' > $CMD_FILE

#setting up the vmd visuals
echo "axes location off" >> $CMD_FILE
echo "color Display Background black" >> $CMD_FILE
echo "mol modstyle 0 0 NewCartoon 0.300000 10.000000 4.100000 0" >> $CMD_FILE
echo "mol modcolor 0 0 Fragment" >> $CMD_FILE
echo "color scale method RGB" >> $CMD_FILE
echo "animate goto 0" >> $CMD_FILE
echo "AutoScaleAllVisible 0.93" >> $CMD_FILE


COUNTER=0

for i in `ls $PDB_FILES`
do

	id=$(printf "%08d" $COUNTER)

	echo "animate goto $COUNTER" >> $CMD_FILE

	echo "render snapshot $TMP_PATH/$id.tga" >> $CMD_FILE

	((COUNTER+=1))
done


echo "quit" >> $CMD_FILE


#run vmd without the GUI, giving it the command file to run
vmd -dispdev openglpbuffer -size $VIDEO_WIDTH $VIDEO_HEIGHT -e $CMD_FILE $PDB_FILES


#file to give to ffmpeg, stating ordering and repeating of frames (1s freezeframe, forwards, 1s freezeframe, backwards)
CONCAT_FILE=$TMP_PATH/concat


echo "file $TMP_PATH/00000000.tga" > $CONCAT_FILE
echo "duration 1" >> $CONCAT_FILE


((COUNTER-=2))

for ((i=1; i<=COUNTER; i++))
do

	id=$(printf "%08d" $i)

	echo "file $TMP_PATH/$id.tga" >> $CONCAT_FILE
	echo "duration $FRAME_LENGTH" >> $CONCAT_FILE

done


echo "file $TMP_PATH/$(printf "%08d" $((COUNTER+1))).tga" >> $CONCAT_FILE
echo "duration 1" >> $CONCAT_FILE


for ((i=COUNTER; i>=0; i--))
do

	id=$(printf "%08d" $i)

	echo "file $TMP_PATH/$id.tga" >> $CONCAT_FILE
	echo "duration $FRAME_LENGTH" >> $CONCAT_FILE

done


#encode video from frames (tgas) using ffmpeg
ffmpeg -f concat -safe 0 -i $TMP_PATH/concat -framerate $FPS -r $FPS -c:v libx264 -pix_fmt yuv420p -crf 20 -maxrate 6M -bufsize 1M -threads 1 -y $DEST_FILE


#remove all the temporary files/folders made in this script, leaving only the video
rm -rf $TMP_PATH