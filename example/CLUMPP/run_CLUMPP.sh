#!/bin/bash

# script to run clump for E. lepidum env, meristic, and morphometric data
# runs clump on all files in directory ending in "paramfile"

# make sure CLUMPP files have UNIX line endings
dos2unix *.paramfile
dos2unit *.indfile

ls *.paramfile | uniq > temp1
FILELIST=`cat temp1`

for f in ${FILELIST}
do
CLUMPP ${f}
done

rm temp1
