#MAPPING OF LIP-MS PEPTIDE SEQUENCE DATA ON PROTEIN STRUCTURES
version 0.1

##How to run the demo:
```
./run_demo.sh input/myoglobin_lip_peptides.tsv input/horse_proteome.fasta output/ /opt/local/bin/pymol
```

or simply

```
./run_demo.sh
```

##Introduction
Limited proteolysis coupled to mass spectrometry (LiP-MS) is a novel proteomics approach that enables the identification of protein structural changes in vivo on a proteome-wide scale. For understanding the functional origin of the LiP patterns, peptides identified in LiP-MS can be mapped on protein structures and their proximity to functionally important sites can be assessed. The mapping will reveal the most likely functional origin of the detected conformational change and allow the distinction of conformotypic peptides from false positive identifications. Molecular viewer like [PyMOL](#pymol) and [VMD](www.ks.uiuc.edu/Research/vmd) can help in interpreting the mapping results. 

The mapping of LiP peptides on protein structures consists of three main steps:

1. Identify the proteins from which the peptides originate using the PERL script [getPeptidePos.pl](#getPeptidePos.pl)
2. Map peptides on protein structures using the PERL script [mapPeptideOnPDB.pl](#mapPeptideOnPDB.pl)
3. Visualize the mapping using your favorite molecular viewer like [PyMOL](#pymol).

	
## Demo Run
This directory has following file structure:

* demo/
   * example_output/
     * *A sample of a finished peptide mapping run*.
   * input/
     * *Input files required for performing the mapping calculations*.
   * output/
     * *The directory in which the results of a new peptide mapping calculation will be stored. Should be empty initially*. 
   * README.md
     * *This read me file*.
   * run\_demo.sh
     * *A master-SHELL-script executing automatically each of the scripts in the scripts directory. Type ```./run_demo.sh``` to run the demo and to get more information*.
   * scripts/
     * *A set of PERL scripts performing each of the steps above*.

The set of scripts performing the peptide mapping on a protein structure are provided in the
scripts directory and can be sequentially executed with the ```./run_demo.sh``` master script. To run the
master script, you will need a working internet connections. If available the path to a local PyMOL 
installation can also be provided to automatically start PyMOL after the mapping has finished. 
The list of the required software tools together with their download link can be found in the table below.
Please type

```
./run_demo.sh
```
to run the demo.

Individual parameters in the master script can be changed within the run_demo.sh script.
The total run time of the demo input data and calculation is around 10 seconds on a modern computer.

The result will be stored in the pml directory within the output directory. The subdirectory should
hold five pml PyMOL scripts for peptide mapping on five crystal structures. The scripts can be loaded
into PyMOL either as an option at the start of PyMOL on the commandline or by using the @-prefix on the 
commandline bar within the PyMOL viewer window.

## System Requirements
- Unix like operating system, recommended and tested on a Mac OS X 10.11.4.
- bash or csh shell
- PERL version 5

##Required Software
<a name="pymol"></a>

Name  | Description                       | Version | URL                                   
------|-----------------------------------|---------|-----
PyMOL | Molecular visualization software. | 1.8     | [http://www.pymol.org](http://www.pymol.org) <br> [http://sourceforge.net/projects/pymol](http://sourceforge.net/projects/pymol)
      | Open-source version of PyMOL can also be installed via MacPorts or RPM | | [https://www.macports.org](https://www.macports.org) <br> [http://rpm.org](http://rpm.org)


##Input files
This demo requires two input files:

1. A tab delimited text file in which the 3rd column holds a list of simple peptide sequences.
2. A standard fasta sequence file of a proteome of interest in which the peptides sequences will be searched. Proteome fasta files can be downloaded from the [UniProt FTP site](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/).

##Commandline options for the scripts
<a name="getPeptidePos.pl"></a>
###getPeptidePos.pl

Option   | required or optional | Default | Information |
---------|----------------------|---------|-------------|
-fasta   | required             |         | Path to a FASTA formatted proteome sequence text file.
-in      | required             |         | Path to a tab delimited text file, which holds in one column a list of peptide sequences.
-col     | optional             | 1       | Column number of peptide sequences in -in file. Numbering starts with 1.
-header  | optional             |         | Ignore header in peptide sequence file.
-verbose | optional             |         | Outputs additional information such as a table header.
-help    | optional             |         | Print this help information.


<a name="mapPeptideOnPDB.pl"></a>
###mapPeptideOnPDB.pl
Option       | required or optional | Default  | Information |
-------------|----------------------|----------|-------------|
-in          | required             |          | Path to a tab delimted text file with peptide and UniProt information as created by  [getPeptidePos.pl](#getPeptidePos.pl)
-del         | optional             | \t       | Column delimiter in -in file.                              
-uniprotCol  | optional             | 1        | Column number of UniProt IDs in input file. Numbering starts with 1.
-pepCol      | optional             | 2        | Column number of peptide sequences in input file. Numbering starts with 1. In case its identical to uniprotCol, split based on _ will be performed and the 2nd string will be taken.
-uniprotDir  | optional             | uniprot  | Directory in which UniProt files will be downloaded. Download of UniProt files already existent in this directory will be skipped.
-fastaDir    | optional             | fasta    | Directory in which protein FASTA files will be written. Download of FASTA files already existent in this directory will be skipped.
-pdbDir	     | optional             | pdb      | Directory in which PDB files will be downloaded. Download of PDB files already existent in this directory will be skipped.
-modelDir    | optional             | homology | Directory in which MODBASE protein model files will be downloaded. Download of MODBASE files already existent in this directory will be skipped.
-pmlDir	     | optional             | pml      | Directory in which PyMOL scripts will be written.
-minSeqId    | optional             | 10       | Minimum sequence identity of homology model to its template structure.
-maxDownload | optional             | 5        | Maximum number of PDB structures to be downloaded.
-header      | optional             |          | Skips the first line of the input file.
-verbose     | optional             |          | Output additional information about the calculation progress on the STDERR channel.
-help        | optional             |          | Print this help information.

