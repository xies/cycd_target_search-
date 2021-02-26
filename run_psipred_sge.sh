#!/bin/bash

filename=$(readlink -f ${1}) # Get absolute path if given relative
#echo "Running psipred on all FASTA files round in ${filename}"
# Expects as argument the directory

for seq_file in ${filename}/*.fasta
do
   /home/xies/Code/psipred/BLAST+/runpsipredplus ${seq_file}
done
