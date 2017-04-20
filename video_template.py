#run this script by typing pymol view.py at the command line in the directory that contains the pdb files you want to view as a movie.
#edit name appropriately - this is for the standard input_froda_* files
from pymol import cmd
#this is to view all of the pdb files in your current working directory in number order, using pymol.
from glob import glob
# print sys.argv
folder=sys.argv[1]
prot=folder.rsplit("/",1)[1]
lst = glob(folder+"/*.pdb")
lst.sort()
for fil in lst: cmd.load(fil,"mov")
cmd.hide()
cmd.show('cartoon')
cmd.bg_color('white')
