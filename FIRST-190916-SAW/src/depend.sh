#!/bin/sh

# get the compiler
CPP="$1"

# get the name of the directory
DIR="$2"
shift 1

# process files in the same directory as the makefile
case "$DIR" in 
"" | "." )
$CPP -MM -MG "$@" | 
sed -e 's@^\(.*\)\.o:@\1.d \1.o:@'
;;
*)
# process files in the module subdirectories. Note the double quotes
# for variable expansion
$CPP -MM -MG "$@" | 
sed -e "s@^\(.*\)\.o:@$DIR/\1.d $DIR/\1.o:@"
;;
esac
