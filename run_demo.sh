#!/bin/bash
if [ $# -lt 3 ]
then
    peptideFile="input/myoglobin_lip_peptides.tsv"
    proteomeFile="input/horse_proteome.fasta"
    outDir="output"
    pymol="/opt/local/bin/pymol"
    peptidePosFile="$outDir/myoglobin_lip_peptides_pos.tsv"

echo
    echo "Demo execution for mapping peptides identified in a LiP-SRM experiment on protein structures. If available, you can provide the path to a PyMOL executable to automatically start PyMOL after the computation has finished."
    echo
    echo "Usage: $0 1-peptideSequenceFile 2-proteomeFastaFile 3-outputDir (4-pymolPath)"
    echo
    echo "Example: $0 input/myoglobin_lip_peptides.tsv input/horse_proteome.fasta output/ /opt/local/bin/pymol"
    echo
    echo "##############################################################################"
    echo "WARNING: No input parameter detected. Running calculations on demo samples:"
    echo "$0 $peptideFile $proteomeFile $outDir $pymol"
    echo 

    while true; do
	read -p "Do you wish to continue ? [yN]: " yn
	case $yn in
            [Yy]* ) echo -e "\nRUNNING DEMO ON DEMO INPUT FILES ...\n"; break;;
            [Nn]* ) echo -e "\nPlease provide costum input parameters as exemplified above.\n"; exit;;
            * ) echo "Please answer yes or no.";;
	esac
    done
else
    peptideFile=$1
    proteomeFile=$2
    outDir=$3
    pymol=$4
fi

# check that all files and directories exist
exit=0
if [ ! -e $peptideFile ]; then echo "    ERROR: Peptide sequence file $peptideFile does not exist. Please verify that it exists."; exit=1; fi 
if [ ! -e $proteomeFile ]; then echo "    ERROR: Proteome fasta file $proteomeFile does not exist. Please verify that it exists."; exit=1; fi 
if [[ ! `type perl` ]]; then echo "    ERROR: PERL could not be found on your operating system. Please download and install it from http://www.perl.org"; exit=1; fi
if [ $exit -eq 1 ]; then exit; fi

# get the absolute paths for all files and directories.
pwd=`pwd`
if [ ! -d $outDir ]; then mkdir -p $outDir; fi
cd $outDir; outDir=`pwd`; cd $pwd

if [ ! -z "$pymol" ] && [ ! -e $pymol ]; then
    echo
    echo "    ########################################################################################################################"
    echo "    WARNING: PyMOL could not be found at $pymol. Please check your path or download and install it from http://www.pymol.org";
    echo "    ########################################################################################################################"
    echo
fi;

echo
echo "2. Identifying the associated proteins of detected peptides"
perl scripts/getPeptidePos.pl -in $peptideFile -fasta $proteomeFile -col 2 -v -header > $peptidePosFile
echo "Done"
echo
echo "1. Mapping peptides on protein structures"
perl scripts/mapPeptideOnPDB.pl -in $peptidePosFile -uniprotCol 2 -pepCol 1 -uniprotDir $outDir/uniprot -fastaDir $outDir/fasta -pdbDir $outDir/pdb -modelDir $outDir/homology -pmlDir $outDir/pml -maxDownload 5 -minSeq 10 -verbose 
echo "-------"
echo "0. DONE"

if [ ! -z "$pymol" ] && [ -e $pymol ]; then
    pmlFile=`find $outDir/pml -name '*.pml' | head -1`
    echo "-------"
    echo "STARTING PyMOL and loading script $pmlFile"
    echo "-------"
    $pymol $pmlFile &
fi
