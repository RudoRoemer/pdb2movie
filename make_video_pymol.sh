#!/bin/bash

#check that five arguments have been given
if [ "$#" -ne 8 ] && [ "$#" -ne 9 ]
then
  echo "usage: $0 video_width video_height frame_rate path_with_pdb_files destination_file_name_and_path video_codec confs freq [pymol_command_file]"
  exit -1
fi

#arguments:
#PDB_PATH is where to find the PDBs which must begin 'tmp_froda_'.
#DEST_FILE is the name (and location) of the video to generate
#OPTIONAL_COMMAND_FILE is optional, and allows extra instructions to be given to pymol before generation of frames
VIDEO_WIDTH=$1
VIDEO_HEIGHT=$2
FPS=$3
PDB_PATH=$4
DEST_FILE=$5
VIDEO_CODEC=$6
CONFS=$7
FREQ=$8
OPTIONAL_COMMAND_FILE=$9

TMP_PATH=$PDB_PATH/pymol_temp

#calculate how many seconds each frame lasts, rounded to ensure ffmpeg does the right thing
FRAME_LENGTH=$(awk -v fps=$FPS 'BEGIN { printf "%.6f", 1.0 / fps }')

#make a temporary folder for the temporary files
rm -rf $TMP_PATH
mkdir $TMP_PATH

#command file to give to pymol, to generate all necessary frames (in png format)
CMD_FILE=$TMP_PATH/pymol.py

PDB_FILES=$PDB_PATH/tmp_froda_*.pdb

#array of all pdb files
arr=(`ls $PDB_FILES`)

###
echo "make_video_pymol.sh: generating command file to give pymol - ${PDB_PATH}"
####

#setting up the pymol visuals
echo 'from pymol import cmd' > $CMD_FILE
echo 'import sys' >> $CMD_FILE

echo "cmd.load(\"${arr[0]}\", 'original')" >> $CMD_FILE
echo 'cmd.spectrum()' >> $CMD_FILE
echo "cmd.viewport($VIDEO_WIDTH, $VIDEO_HEIGHT)" >> $CMD_FILE
echo 'cmd.hide()' >> $CMD_FILE
echo 'cmd.show("cartoon")' >> $CMD_FILE
echo 'cmd.bg_color("grey")' >> $CMD_FILE

#if an optional file with extra instructions for pymol has been included, use it
echo "OPTIONAL_COMMAND_FILE:" $OPTIONAL_COMMAND_FILE
if [ $OPTIONAL_COMMAND_FILE != "" ]
then
  cat $OPTIONAL_COMMAND_FILE >>  $CMD_FILE
fi

echo "cmd.save(\"$PDB_PATH/cartoon.pse\")" >> $CMD_FILE

COUNTER=0

for i in `ls $PDB_FILES`
do

  #get the filename, and extract the conformer number
  j="$(basename -- $i)"
  j=${j:10:8}
  if [ "$j" -eq "00000000" ]
  then
    j=0
  else
    j=$(echo $j | sed 's/^0*//')
  fi

  #check if we want to use this pdb; if not, then skip it

  #its number must be less than CONFS...
  if (( $j > $CONFS ))
  then
    continue
  fi

  #...and must be a multiple of FREQ
  if (( $j % $FREQ != 0 ))
  then
    continue
  fi

  #id of the next frame (png file) to generate
  id=$(printf "%08d" $COUNTER)

  #copy the original (initial) pdb to mov
  echo 'cmd.copy("mov", "original")' >> $CMD_FILE

  #load the next pdb into mov
  echo "cmd.load(\"$i\", 'mov')" >> $CMD_FILE

  #display every currently loaded object in cartoon form
  echo 'cmd.show("cartoon")' >> $CMD_FILE

  #hide (undisplay) the original (initial) pdb
  echo 'cmd.hide("everything", "original")' >> $CMD_FILE

  #go to frame 2 of mov (the frame containing the pdb to be rendered next)
  echo "cmd.frame(2)" >> $CMD_FILE

  #render the next frame (png file)
  echo "cmd.mpng(\"$TMP_PATH/$id\", 2, 2, 0, 0, 2, 0, $VIDEO_WIDTH, $VIDEO_HEIGHT)" >> $CMD_FILE

  #remove (unload) mov (including the pdb just used)
  echo 'cmd.delete("mov")' >> $CMD_FILE

  ((COUNTER+=1))
done

echo 'cmd.quit()' >> $CMD_FILE

###
echo "make_video_pymol.sh: running pymol with generated command file to create frames (png files) - ${PDB_PATH}"
echo "make_video_pymol.sh: copy of pymol command file following NOW"
cat $CMD_FILE
###

#run pymol without the GUI, suppressing the startup output, giving it the command file to run
( set -x ; pymol -c -q $CMD_FILE )

###
echo "make_video_pymol.sh: generating concat file to give ffmpeg - ${PDB_PATH}"
###

#file to give to ffmpeg, stating ordering and repeating of frames (1s freezeframe, forwards, 1s freezeframe, backwards)
CONCAT_FILE=$TMP_PATH/concat


#hold the first png file for one second
echo "file $TMP_PATH/000000000002.png" > $CONCAT_FILE
echo "duration 1" >> $CONCAT_FILE

#use each png file in forward order, holding each for the length of a frame
((COUNTER-=2))

for ((i=1; i<=COUNTER; i++))
do

  id=$(printf "%08d" $i)

  echo "file $TMP_PATH/${id}0002.png" >> $CONCAT_FILE
  echo "duration $FRAME_LENGTH" >> $CONCAT_FILE

done

#hold the last png file for one second
echo "file $TMP_PATH/$(printf "%08d" $((COUNTER+1)))0002.png" >> $CONCAT_FILE
echo "duration 1" >> $CONCAT_FILE

#use each png file in reverse order, holding each for the length of a frame
for ((i=COUNTER; i>=0; i--))
do

  id=$(printf "%08d" $i)

  echo "file $TMP_PATH/${id}0002.png" >> $CONCAT_FILE
  echo "duration $FRAME_LENGTH" >> $CONCAT_FILE

done

###
echo "make_video_pymol.sh: running ffmpeg with generated concat file to create video from frames using $VIDEO_CODEC codec - ${PDB_PATH}"
###

#encode video from frames (png files) using ffmpeg and the specified codec

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

#remove all the temporary files/folders made in this script, leaving only the video
rm -rf $TMP_PATH
