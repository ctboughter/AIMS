#!/bin/bash

# So this script takes *csv files and makes them into readable FASTA formats for igblast
# Should also extract relevant biophysical properties for these sequences while we're at it

# ASSUME FOR NOW THAT EVERY XLSX FILE HAS A HEADER. Should (maybe?) be some way to detect this
# Honestly detection should be hard, should probably put some kind of IF statement in...

# If we're using CSV, this is the general way to pull out columns
#awk -F ',' '{print $1}' drug_seqs.csv

# Lastly, make sure that you're running these scripts for the first time...

# Alright let's get into it
echo 'Which file has the amino acid sequences? (full name, like test.txt)'
read seqs_file

echo 'In what column can you find the heavy chain sequences? (1,2,3...)'
read heavy_col

echo 'In what column can you find the light chain sequences? (1,2,3...)'
read light_col

awk -v tt=$heavy_col -F ',' '{print $tt}' $seqs_file > testH.txt

for i in $(cat testH.txt); do echo '>' >> seqs_fastaH.txt; echo $i >> seqs_fastaH.txt; done
for i in $(cat testH.txt); do echo $i >> seqsH.txt; done

awk -v tt=$light_col -F ',' '{print $tt}' $seqs_file > testL.txt

for i in $(cat testL.txt); do echo '>' >> seqs_fastaL.txt; echo $i >> seqs_fastaL.txt; done
for i in $(cat testL.txt); do echo $i >> seqsL.txt; done

rm testH.txt testL.txt

# Should add in an extra line to remove the header, but can probably also do this automatically...