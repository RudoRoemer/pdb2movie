#!/bin/bash

#check that five arguments have been given
if [ "$#" -ne 6 ] && [ "$#" -ne 7 ]
then
  echo "usage: $0 video_width video_height frame_rate path_with_pdb_files destination_file_name_and_path video_codec [pymol_command_file]"
  exit -1
fi

#arguments:
#PDB_PATH is where to find the PDBs which must begin 'tmp_froda_'.
#DEST_FILE is the name (and location) of the video to generate
#OTIONAL_COMMAND_FILE is optional, and allows extra instructions to be given to pymol before generation of frames
VIDEO_WIDTH=$1
VIDEO_HEIGHT=$2
FPS=$3
PDB_PATH=$4
DEST_FILE=$5
VIDEO_CODEC=$6
OTIONAL_COMMAND_FILE=$6

TMP_PATH=$PDB_PATH/vmd_temp

#calculate how many seconds each frame lasts, rounded to ensure ffmpeg does the right thing
FRAME_LENGTH=$(awk -v fps=$FPS 'BEGIN { printf "%.6f", 1.0 / fps }')

#make a temporary folder for the temporary files
rm -rf $TMP_PATH
mkdir $TMP_PATH

#command file to give to vmd, to generate all necessary frames (in tga format)
CMD_FILE=$TMP_PATH/vmd.cmd

PDB_FILES=$PDB_PATH/tmp_froda_*.pdb

###
echo "make_video_vmd.sh: generating command file to give vmd - ${PDB_PATH}"
####

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
echo "animate goto 0" >> $CMD_FILE
echo "axes location off" >> $CMD_FILE
echo "color Display Background black" >> $CMD_FILE
echo "mol modstyle 0 0 NewCartoon 0.300000 10.000000 4.100000 0" >> $CMD_FILE
echo "mol modcolor 0 0 Fragment" >> $CMD_FILE
echo "color scale method RGB" >> $CMD_FILE
echo "AutoScaleAllVisible 0.93" >> $CMD_FILE

#if an optional file with extra instructions for vmd has been included, use it
if [ "$OPT_SETTINGS_FILE" != "" ]
then
  cat $OTIONAL_COMMAND_FILE >>  $CMD_FILE
fi

COUNTER=0

for i in `ls $PDB_FILES`
do
  #id of the next frame (tga file) to generate
  id=$(printf "%08d" $COUNTER)

  #load the next pdb
  echo "mol addfile {$i} type {pdb} first 0 last -1 step 1 waitfor 1 0" >> $CMD_FILE

  #goto the pdb just loaded
  echo "animate goto 1" >> $CMD_FILE

  #render a tga file (the next frame)
  echo "render snapshot $TMP_PATH/$id.tga" >> $CMD_FILE

  #remove (deload) the pdb just used
  echo "animate delete beg 1 end 1 skip 0 0" >> $CMD_FILE

  ((COUNTER+=1))
done

echo "quit" >> $CMD_FILE


#array of all pdb files
arr=(`ls $PDB_FILES`)

###
echo "make_video_vmd.sh: running vmd with generated command file to create frames (tga files) - ${PDB_PATH}"
###

#run vmd without the GUI, giving it the command file to run
( set -x ; vmd -dispdev openglpbuffer -size $VIDEO_WIDTH $VIDEO_HEIGHT -e $CMD_FILE ${arr[0]} )

###
echo "make_video_vmd.sh: generating concat file to give ffmpeg - ${PDB_PATH}"
###

#file to give to ffmpeg, stating ordering and repeating of frames (1s freezeframe, forwards, 1s freezeframe, backwards)
CONCAT_FILE=$TMP_PATH/concat


#hold the first png file for one second
echo "file $TMP_PATH/00000000.tga" > $CONCAT_FILE
echo "duration 1" >> $CONCAT_FILE

#use each png file in forward order, holding each for the length of a frame
((COUNTER-=2))

for ((i=1; i<=COUNTER; i++))
do

  id=$(printf "%08d" $i)

  echo "file $TMP_PATH/$id.tga" >> $CONCAT_FILE
  echo "duration $FRAME_LENGTH" >> $CONCAT_FILE

done

#hold the last png file for one second
echo "file $TMP_PATH/$(printf "%08d" $((COUNTER+1))).tga" >> $CONCAT_FILE
echo "duration 1" >> $CONCAT_FILE

#use each png file in reverse order, holding each for the length of a frame
for ((i=COUNTER; i>=0; i--))
do

  id=$(printf "%08d" $i)

  echo "file $TMP_PATH/$id.tga" >> $CONCAT_FILE
  echo "duration $FRAME_LENGTH" >> $CONCAT_FILE

done

###
echo "make_video_vmd.sh: running ffmpeg with generated concat file to create video from frames - ${PDB_PATH}"
###

#encode video from frames (tga files) using ffmpeg and the specified codec

if [ "$VIDEO_CODEC" == "mp4" ]
then
  (
  set -x
  
  ffmpeg -hide_banner -loglevel warning -f concat -safe 0 -i $TMP_PATH/concat -framerate $FPS -r $FPS -c:v libx264 -pix_fmt yuv420p -crf 20 -maxrate 6M -bufsize 1M -threads 1 -movflags faststart -y ${DEST_FILE}
  )
fi

if [ "$VIDEO_CODEC" == "hevc" ]
then
  (
  set -x
  
  ffmpeg -hide_banner -loglevel warning -f concat -safe 0 -i $TMP_PATH/concat -framerate $FPS -r $FPS -c:v libx265 -pix_fmt yuv420p -crf 25 -maxrate 6M -bufsize 1M -threads 1 -movflags faststart -y ${DEST_FILE}
  )
fi


#set the correct permissions
chmod 744 $DEST_FILE


#remove all the temporary files/folders made in this script, leaving only the video
rm -rf $TMP_PATH