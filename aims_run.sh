#!/bin/sh

# This file solely exists so that you can keep track of the MANY options for the AIMS CLI
# Define which example analysis you want to run OR comment out the other 

##############################################################
# Define the molecule as either ab, tcr, peptide, or msa
# FOR BASH, CANT HAVE A SPACE BETWEEN "=" signs
molecule=msa
##############################################################

# An example of how to run the CLI for TCR analysis
# Only show one file here as an example, but it works with as many
# files as you want to look at.
if [ $molecule = ab ]
then
    python aims_cli.py \
    --datDir app_data \
    --outputDir AIMS_ab \
    --fileNames flu_mono.csv \
    --datNames flu \
    --numLoop 6 \
    --molecule ig > aims_ab.out
    #Alternative option 
    #--fileNames siv_tl8.csv siv_cm9.csv \
    #--datNames #TL8 CM9 \
fi

if [ $molecule = tcr ]
then
    python aims_cli.py \
    --datDir app_data \
    --outputDir AIMS_tcr \
    --fileNames siv_tl8.csv siv_cm9.csv \
    --datNames TL8 CM9 \
    --numLoop 1 \
    --molecule ig > aims_tcr.out
fi

if [ $molecule = peptide ]
then
    python aims_cli.py \
    --molecule peptide \
    --datDir app_data \
    --outputDir AIMS_pep \
    --fileNames pancreas_hla_atlas.csv kidney_hla_atlas.csv \
    --datNames Pancreas Kidney > aims_pep.out
fi
# Running the CLI for MSA/MHC Analysis
# Because the MSA analysis is a bit complicated, go all out and use this as an opportunity 
# to show how almost ALL of the CLI options work. Users could then copy and paste these
# into their own custom CLI commands as they see fit.

if [ $molecule = msa ]
then
    python aims_cli.py \
    --datDir app_data \
    --outputDir AIMS_mhc \
    --molecule MSA \
    --fileNames cd1.fasta classIa.fasta fish.fasta \
    --datNames CD1 ClassIa Fish \
    --align center \
    --subset True \
    --subStart 164 214 275 327 \
    --subEnd 214 275 327 376 \
    --dropDup True \
    --normProp True \
    --clustData avg \
    --projAlg pca \
    --umapSeed 42 \
    --clustAlg kmean \
    --clustSize 3 \
    --metaForm category \
    --metaName Dset \
    --showProj both \
    --showClust both \
    --normBar False \
    --analysisSel metadata \
    --saveSeqs True \
    --selDat 0 1 2 \
    --prop1 3 \
    --prop2 4 \
    --matSize 5 > aims_msa.out
fi

