#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 26.08.2016

use strict;
use warnings;
use Cwd;
use File::Basename;
use Getopt::Long;
use File::Path;
use File::Fetch;
use List::Util qw[min max sum];
use Scalar::Util qw(looks_like_number);
use File::Copy;
use File::Path qw(make_path);

# do not buffer output
$| = 1;

##############################################################################
### Calculates various properties for a peptide such as solvent            ###
### temperature factor, occupancy value, gap size, distance to 2nd         ###
### structure elements etc.                                                ###
##############################################################################
my (
    # variable for parameters which are read in from commandline
    $help,

    $infile,
    $del,
    $uniprotCol,
    $peptideCol,

    $uniprotDir,
    $fastaDir,
    $pdbDir,
    $modelDir,
    $pmlDir,

    $minModBaseSeqId,
    $maxDownload,

    $verbose,
    $hasHeader,
    $xDelTrash,

   );

#################
# global variable
my $nontrypticSelection = "";
my $nontrypticName = "nontryptic_is_red";
##############################################################################
### read all needed parameters from commandline ##############################
##############################################################################
my @argv = @ARGV;

&GetOptions(
    "help!"          => \$help,            # Print this help.
    "in=s"           => \$infile,          # Tab delimted text file with peptide and UniProt information.
    "del:s"          => \$del,             # Delimiter if -in file is not tab delimited.
    "uniprotCol:s"   => \$uniprotCol,      # Column number of UniProt IDs e.g. KAR9_YEAST, in input file. Column numbering starts with 1. Default is 2.
    "pepCol:s"       => \$peptideCol,      # Column number of peptide sequence in input file. Column numbering starts with 1. Default is 1. In case its identical to uniprotCol, split based on _ will be performed and the 2nd string will be taken  
    "uniprotDir:s"   => \$uniprotDir,      # Directory in which UniProt files will be written. Download of UniProt files already existent in this directory will be skipped. Default is ./uniprot 
    "fastaDir:s"     => \$fastaDir,        # Directory in which protein FASTA files will be written. Download of FASTA files already existent in this directory will be skipped. Default is ./fasta 
    "pdbDir:s"       => \$pdbDir,          # Directory in which PDB files will be written. Download of PDB files already existent in this directory will be skipped. Default is ./pdb
    "modelDir:s"     => \$modelDir,        # Directory in which MODBASE protein model files will be written. Download of MODBASE files already existent in this directory will be skipped. Default is ./homology
    "pmlDir:s"       => \$pmlDir,          # Directory in which PyMOL scripts will be written. Only valid in combination with -pymol parameter. Default is ./pml
    "minSeqId:i"     => \$minModBaseSeqId, # Minimum sequence identity of homology model to its template structure, in order to be considered for the peptide mapping. Default is 10.
    "maxDownload:i"  => \$maxDownload,     # Maximum number of PDB structures to be downloaded for mapping. Default is 5.
    "header!"        => \$hasHeader,       # Skips the first line of the input file.
    "verbose!"       => \$verbose,         # Print out additional information about the calculation progress on the STDERR channel.
) or die "\nTry \"$0 -h\" for a complete list of options\n\n";

##############################################################################
# HELP
if ($help) {printHelp(); exit}
##############################################################################
# INITIAL GLOBAL SETTINGS

$del = "\t" if(!defined $del);
$maxDownload = 5 if(!defined $maxDownload);

$minModBaseSeqId = 10 if(!defined $minModBaseSeqId);

my $pwd = getcwd;
$pwd = &fullPath($pwd);
my $uniprotSeqAlignFile = "uniprot";
my $pdbSeqAlignFile = "pdb";

my $fileEndUniProt = ".txt";
my $fileEndFasta = ".fasta";
my $fileEndPdb = ".pdb1";

$uniprotCol = 2 if(!defined $uniprotCol);
$uniprotCol--;
$peptideCol = 1 if(!defined $peptideCol);
$peptideCol--;

$uniprotDir  = "$pwd/uniprot"   if(!defined $uniprotDir);
$fastaDir    = "$pwd/fasta"     if(!defined $fastaDir);
$pdbDir      = "$pwd/pdb"       if(!defined $pdbDir);
$modelDir    = "$pwd/homology"  if(!defined $modelDir);
$pmlDir      = "$pwd/pml"       if(!defined $pmlDir);

my %aa1 = ("G" => "GLY",
	   "A" => "ALA", 
	   "V" => "VAL",
	   "L" => "LEU",
	   "I" => "ILE",
	   "F" => "PHE",
	   "P" => "PRO",
	   "W" => "TRP",
	   "N" => "ASN",
	   "Q" => "GLN",
	   "M" => "MET",
	   "C" => "CYS",
	   "T" => "THR",
	   "Y" => "TYR",
	   "S" => "SER",
	   "R" => "ARG",
	   "K" => "LYS",
	   "H" => "HIS",
	   "D" => "ASP",
	   "E" => "GLU",
	   "B" => "ASX",
	   "Z" => "GLX");

my %aa3 = ("GLY" => "G",
	   "ALA" => "A", 
	   "VAL" => "V",
	   "LEU" => "L",
	   "ILE" => "I",
	   "PHE" => "F",
	   "PRO" => "P",
	   "TRP" => "W",
	   "ASN" => "N",
	   "GLN" => "Q",
	   "MET" => "M",
	   "CYS" => "C",
	   "THR" => "T",
	   "TYR" => "Y",
	   "SER" => "S",
	   "ARG" => "R",
	   "LYS" => "K",
	   "HIS" => "H",
	   "ASP" => "D",
	   "GLU" => "E",
	   "ASX" => "B",
	   "GLX" => "Z");

my $blosum62 = 
"#AA A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X -\n".
"A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4\n".
"R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4\n".
"N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4\n".
"D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4\n".
"C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4\n".
"Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4\n".
"E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n".
"G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4\n".
"H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4\n".
"I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4\n".
"L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4\n".
"K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4\n".
"M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4\n".
"F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4\n".
"P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4\n".
"S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4\n".
"T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4\n".
"W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4\n".
"Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4\n".
"V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4\n".
"B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4\n".
"Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n".
"X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4\n".
"- -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1\n";
my %h = ();
my %blosum62 = ();
my @l = split(/\n/, $blosum62);
for(my $n = 0; $n < @l; $n++) {
    my @a = split(/\s+/, $l[$n]);
    for(my $i=0; $i < @a; $i++) {
	if($l[$n] =~ /^#/) {
	    $h{$i} = $a[$i];
	}
	else{
	    $blosum62{$h{$n}}->{$h{$i}} = $a[$i];
	}
    }
}


##############################################################################
##############################################################################
### SUBROUTINES
##############################################################################

##############################################################################
###############################################################################
sub printHelp {
###############################################################################
    # prints a help about the using and parameters of this scripts 
    # (execute if user types commandline parameter -h)
    # param:  no paramaters
    # return: no return value

    my (
	$usage,
	$sourceCode,
	@rows,
	$row,
	$option,
	$scriptInfo,
	$example,
       );

    $usage = "$0 -in myoglobin_lip_peptides_pos.tsv -uniprotCol 2 -pepCol 1 -header -v";

    print "\nUsage: " .  $usage . "\n\n";

    print "Valid options are:\n\n";
    open(MYSELF, "$0") or
      die "ERROR: Cannot read source code file $0: $!.\n";
    $sourceCode .= join "", <MYSELF>;
    close MYSELF;
    $sourceCode =~ s/^.+?\&GetOptions\(\n//s;
    $sourceCode =~ s/\n\).+$//s;
    @rows = split /\n/, $sourceCode;
    foreach $row (@rows){
        $option = $row;
	$option =~ s/\s+\"//g;
	$option =~ s/\"\s.+\#/\t\#/g;
	$option =~ s/=./\t<value> [required]/;
	$option =~ s/:./\t<value> [optional]/;
	$option =~ s/!/\t<non value> [optional]/;

	$row =~ s/^.*//;
	print "\t";
	printf("%-1s%-30s%-30s\n", "-",$option,$row);

    } # end of foreach $row (@rows)
    print "\n";
    print "Options may be abreviated, e.g. -h for --help\n\n";

    $example  = "$0";
}
################################################################################
sub checkInputParameter {
    
    if($#argv == -1){
	&printHelp();
	exit;
    } 

    my $doExit = 0;
    if(!defined $infile){
	print STDERR "\nERROR: Please specify the input peptide file. See -help\n";
	$doExit = 1;
    } else {
	if(!-e $infile){
	    print STDERR "\nERROR: Peptide input file could not be found at $infile: $!\n";
	    $doExit = 1;
	} else {
	    # check that columns in input file hold correct information
	    open(F,$infile) or die "\nERROR: Failed to open $infile: $!.\n";
	    my $i = 0;
	    while(my $line = <F>){
		next if($i==0 and $hasHeader);
		$i++;
		# read only first 10 lines
		last if($i==10);
		next if($line=~/^\s*#/);
		chomp($line);
		my @a = split(/$del/,$line);

		my $uniProtId = "";
		my $peptide = "";

		if($uniprotCol == $peptideCol) {
		    ($uniProtId, $peptide, my @tmp) = split(/\_/, uc($a[$uniprotCol])); 
		} else {
		    $uniProtId = uc($a[$uniprotCol]);
		    $peptide = uc($a[$peptideCol]);
		}
		$uniProtId =~ s/.*\|//;

		if($uniProtId !~ /[A-Z0-9]{1,5}_[A-Z0-9]{1,5}/ and $uniProtId !~ /[A-Z][0-9][A-Z0-9]{3}[0-9]/){
		    print STDERR "\nWARNING: $uniProtId in the ".($uniprotCol+1)."th column, line $i does not correspond to an official ".
			         "UniProt ID or Accession number. Did you set the correct column number? See -help\n";
		    next;
		}
		$peptide =~ s/\[.*?\]//g;
		if($peptide =~ /[^A-Z]/){
		    print STDERR "\nWARNING: $peptide in the ".($peptideCol+1)."th column, line $i, seems not to be an amino acid sequence. ".
			         "Did you set the correct column number? See -help\n";
		    next;
		}
	    }
	    close(F);
	}
    }

    exit if($doExit == 1);
}

################################################################################
sub readInputFile {
    my $peptides;
    open(F,$infile) or die "ERROR: Failed to open $infile: $!.\n";
    my $i = 0;
    while(my $line = <F>){
	if(defined $hasHeader) {
	    undef($hasHeader);
	    next;
	}
	next if($line=~/^\s*#/);
	chomp($line);
	my @a = split(/$del/,$line);
	my $uniProtId = "";
	my $peptide = "";

	if($uniprotCol == $peptideCol) {
	    ($uniProtId, $peptide, my @tmp) = split(/\_/, uc($a[$uniprotCol]));
	} else {
	    $uniProtId = uc($a[$uniprotCol]);
	    $peptide = uc($a[$peptideCol]);
	}

	$uniProtId =~ s/.*\|//;
 
	if($uniProtId !~ /[A-Z0-9]{1,5}_[A-Z0-9]{1,5}/ and $uniProtId !~ /[A-Z][0-9][A-Z0-9]{3}[0-9]/){
	    next;
	}
	$peptide =~ s/\[.*?\]//g;
	if($peptide =~ /[^A-Z]/){
	    next;
	}

	$peptides->{$uniProtId}->{$peptide} = ++$i;
    }
    close(F);
    return $peptides;
}
################################################################################
sub changeDir {
    (my $dir) = @_; 
    print STDERR "Changing to directory $dir\n" if(defined $verbose);
    chdir($dir);
}
##############################################################################
sub download {
    (my $url, my $outdir) = @_;

    print STDERR "DOWNLOADING $url ...\n" if(defined $verbose);
    my $ff = File::Fetch->new(uri => $url);
    my $file = $ff->fetch(to => $outdir) or
	print STDERR "ERROR: Failed to download $url: ".$ff->error()."\n";
    return $file;
}


##############################################################################
sub downloadUniProt {
    (my $id, my $fileSuffix, my $outdir) = @_;

    if(-f "$outdir/$id$fileSuffix") {
	    print STDERR "Skipping download as UniProt file already exists at $outdir/$id$fileSuffix\n"
		if(defined $verbose);
	    return "$outdir/$id$fileSuffix";
    }

    if(!-d $outdir){
	print STDERR "Creating directory $outdir\n" if(defined $verbose);
	make_path($outdir);
    }

    my $pageURL = 'http://www.uniprot.org/uniprot/'.$id.$fileSuffix;
    my $file = &download($pageURL, $outdir);
    sleep 1;
    
    if($file) {
	open(F, $file) or die "ERROR: Failed to open downloaded file $file: $!\n";
	my $newId = "";
	my $header = <F>;
	my $isHtml = 1 if($header =~ /DOCTYPE html/);
	if($isHtml) {
	    while(my $line = <F>) {
		if($line =~ /.*\/uniprot\/(.*?)\?.*/){
		    $newId = $1;
		    last;
		}
	    }
	    &downloadUniProt($newId, $fileSuffix, $outdir);
	    print STDERR "Renaming file $newId$fileSuffix to $id$fileSuffix.\n"
		if(defined $verbose);
	    rename($file, "$outdir/$id$fileSuffix");
	    return;
	}
	close(F);
    } else {
	print STDERR "ERROR while downloading UniProt entry from page $pageURL$id$fileSuffix\n";
    }
    return $file;
}
##############################################################################
sub downloadPDB {

    (my $id, my $outdir, my $chain, my $force) = @_;

    if($chain eq "0" or $chain eq "1"){
	if(-f "$outdir/$id$fileEndPdb" and !$force) {
	    print STDERR "Skipping download as PDB file already exists at $outdir/$id$fileEndPdb\n"
		if(defined $verbose);
	    return;
	}
    } else {
	if(-f "$outdir/$id$chain$fileEndPdb" and !$force) {
	    print STDERR "Skipping download as PDB file already exists at $outdir/$id$chain$fileEndPdb\n"
		if(defined $verbose);
	    return;
	}
    }

    if(!-d $outdir){
	print STDERR "Creating directory $outdir\n" if(defined $verbose);
	make_path($outdir);
    }

    my $pageURL = 'http://www.rcsb.org/pdb/files/'.$id.".pdb1.gz";
    my $file = &download($pageURL, $outdir);
    if($file) {
	print STDERR "Unzipping PDB file ... and content to file \"$id$fileEndPdb\"\n"
	    if(defined $verbose);
	`gunzip -f $file`;

	open(F, "$outdir/$id$fileEndPdb") or 
	    die "ERROR: Failed to open \"$outdir/$id$fileEndPdb\": $!.\n";
	my %output = ();
	while(my $l=<F>){
	    if($l =~ /^ENDMDL/){
		last;
	    }
	    if($l =~ /^ATOM/){
		next if(!exists $aa3{substr($l,17,3)});

		substr($l, 21, 1, "A") if(substr($l, 21, 1) eq " ");		
		my $c = substr($l, 21, 1);

		next if($chain ne "0" and $chain ne "1" and 
			$chain =~ /[A-Za-z]/ and $chain !~ /$c/);

		$output{$c} .= $l;
	    } elsif($l =~ /^HETATM/){
		substr($l, 21, 1, "A") if(substr($l, 21, 1) eq " ");		
		$output{HETATM} .= $l;
	    }
	}
	close(F);

	if($chain eq "0" or $chain eq "1"){
	    open(O, ">$outdir/$id$fileEndPdb") or 
		die "ERROR: Failed to create PDB file \"$outdir/$id$fileEndPdb\": $!.\n";
	} else {
	    open(O, ">$outdir/$id$chain$fileEndPdb") or 
		die "ERROR: Failed to create PDB file \"$outdir/$id$chain$fileEndPdb\": $!.\n";
	    print STDERR "Selecting only $chain chain IDs from PDB $id.\n" if(defined $verbose);
	}
	foreach my $chainId (sort {$a cmp $b} keys %output){
           print O $output{$chainId};
            print O "TER\n";
            if($chain eq "1"){
                print STDERR "Removing all chains from \"$id$fileEndPdb\" except of ".
		              "first chain \"$chainId\".\n" if(defined $verbose);
                last;
	    }
	}
	print O $output{HETATM} if(defined $output{HETATM});
	close(O);
    } else {
	print STDERR "ERROR while downloading PDB model from page $pageURL$id.pdb.gz\n";
	return -1;
    }
    return 1;
}
##############################################################################
sub alignmentModBase {
    (my $id, my $outdir) = @_;
    
    return if(-f "$outdir/modbase_"."$id.xml");

    if(!-d $outdir){
	print STDERR "Creating directory $outdir\n" if(defined $verbose);
	make_path($outdir);
    }

    my $pageURL = 'http://salilab.org/modbase/retrieve/modbase/?type=alignment&databaseID='.$id;
    my $file = &download($pageURL, $outdir);

    my $align;
    if($file) {
	print STDERR "Download was successfull.\n"
	    if(defined $verbose);
	open(O, ">$outdir/modbase_"."$id.xml") or die "Couldn't create file \"$outdir/modbase_"."$id.xml\": $!\n";
	open(F, $file) or die "ERROR: Failed to open $file: $!\n";
	my $header;
	my $read = 0;
	while(my $l = <F>){
	    if($l =~ /^>/){
		$read = 1;
		$header = $l;
	    } elsif($l =~ /^structure/ or $l =~ /^sequence/){
		$header .= " $l";
	    } elsif($read) {
		$align->{structure} .= $l if($header=~/structure/);
		$align->{sequence} .= $l if($header=~/sequence/);
		$read = 0 if($l =~ /\*$/);
	    }
	    last if($l =~ /\<\/content\>/);
	    print O "$l\n" if($read);
	}
	close(O);
	close(F);
	unlink $file;
    } else {
	print STDERR "ERROR while downloading MODELLER model from page $pageURL$id\n";
    }
    return $align;
}

##############################################################################
# downloads only the first model from MODBASE to a protein structure.
sub downloadModBase {
    (my $id, my $outdir) = @_;
    
    if(-f "$outdir/modbase_"."$id$fileEndPdb") {
	    print STDERR "SKIPPING DOWNLOAD as MODBASE homology file already exists at $outdir/modbase_"."$id$fileEndPdb\n"
		if(defined $verbose);
	return 1;
    }

    if(!-d $outdir){
	print STDERR "Creating directory $outdir\n" if(defined $verbose);
	make_path($outdir);
    }

    my $align = &alignmentModBase($id, $outdir);
    if((keys %$align) != 2){
	print STDERR "WARNING: Could not find MODBASE alignment for protein $id.\n"
	    if(defined $verbose);
	return 0;
    }
    
    my $pageURL = 'http://salilab.org/modbase/retrieve/modbase/?databaseID='.$id;
    my $file = &download($pageURL, $outdir);

    my $lowQuality = 0;
    if($file) {
	print STDERR "Download was successfull.\n"
	    if(defined $verbose);
	my %uniqAA;
	my $i;
	my $isGap;
	open(O, ">$outdir/modbase_"."$id$fileEndPdb") or die "ERROR: Couldn't create file \"modbase_"."$id$fileEndPdb\": $!\n";
	open(F, $file) or die "ERROR: Couldn't read $file: $!\n";
	while( my $l = <F>){
	    next if($l =~ /\</);
	    last if($l =~ /^END/);
	    if(defined $align->{structure} and defined $align->{sequence}){
		if($l =~ /^REMARK 220 SEQUENCE IDENTITY/){
		    my $seqId = $l;
		    $seqId =~ s/.*://;
		    $seqId =~ s/\s//g;
		    if($seqId < $minModBaseSeqId){
			$lowQuality = $seqId;
			last;
		    }
		}
		if($l =~ /^ATOM  /){
		    my $resName = substr($l, 17, 3);
		    $resName =~ s/\s//g;
		    my $resNo   = substr($l, 22, 4);
		    $resNo =~ s/\s//g;
		    my $chainId = substr($l, 21, 1);

		    my $aa = "$resName-$resNo-$chainId";
		    if(!exists $uniqAA{$aa}){
			$uniqAA{$aa} = $aa;
			$i = keys %uniqAA;
			my $iseq = substr($align->{sequence}, $i-1, 1);
			my $istr = substr($align->{structure}, $i-1, 1);

			$isGap = $iseq eq "-" or $istr eq "-" ? 1 : 0;

#			if($isGap == 0 and $iseq ne $aa3{$resName}) {
#			    print STDERR "WARNING: MODBASE alignment does ".
#				"not match MODBASE PDB sequence: $i: $iseq vs ".
#				"$aa3{$resName} ... skipping residue.\n";
#			}

		    }
		    substr($l, 21, 1, "A") if($chainId eq " ");
		    next if($isGap);
		}
	    }
	    print O "$l\n";
	}
	close(O);
	close(F);

	if($lowQuality > 0){
	    print STDERR "Deleting MODBASE model for $id as its sequence ID to the template is below $minModBaseSeqId"."% with $lowQuality\%.\n"
		if(defined $verbose);
	    unlink("$outdir/modbase_"."$id$fileEndPdb");
	    unlink("$outdir/modbase_"."$id.xml");
	    return 0;
	}

	return 1;
    } else {
	print STDERR "ERROR while downloading MODBASE model from page $pageURL$id\n";
	return 0;
    }
}

################################################################################
sub downloadStructure {
    (my $pdbId, my $outdir) = @_;
    &downloadPDB($pdbId, $outdir, 1, 0);
    return "PDB";
}
################################################################################
sub pdbIdsFromUniProt { 
   (my $uniProtFile) = @_;

   return unless (-e $uniProtFile);

    my $pdbIds;
    my $pdbIdsHighRes;
    my $highResI = 0;
    open(F,$uniProtFile) or die "ERROR: Failed to open $uniProtFile: $!.\n";
    
    while(my $line = <F>){
	chomp($line);
	next if($line=~/^\s*#/);

	if($line =~ /^DR\s+(PDB;.*)/){
	    (my $db, my $pdbId, my $meth, my $res, my $len) = split(/;\s+/, $line);

	    $pdbId = lc($pdbId);

	    # exclude all non X-ray/NMR structures.
	    next if($meth ne "X-ray" and $meth ne "NMR");
	    $res =~ s/ A$//;
	    my $chainIds = $len;
	    $chainIds =~ s/\=.*//;
	    $len =~ s/.*=(.*)\./$1/; 	
	    
	    my @chainIds = split(/\//, $chainIds);
	    foreach my $chainId (@chainIds){
		if(looks_like_number($res) && $res <= 3.0){
		    $pdbIdsHighRes->{$res}->{$pdbId}->{$chainId} = 1;
		    $highResI++;
		}
		$pdbIds->{$res}->{$pdbId}->{$chainId} = 1;
		# use only 1st chain
		last;
	    }
	}
    }
    close(F);
    
    return $highResI > 0 ? $pdbIdsHighRes : $pdbIds;
}
################################################################################
sub pdb2fasta {
    (my $pdbFile) = @_;

    return unless (-e $pdbFile);

    my $pdb = $pdbFile;
    $pdb =~ s/.*\///;
    $pdb =~ s/\..*//;
    my $fasta = "";
    my %uniqAA;
    my $n=-1;
    open(F, $pdbFile) or
	die "ERROR: Failed to open PDB file $pdbFile: $!.\n";
    while(my $l = <F>){
	if($l =~ /^ATOM/){
	    my $resName = substr($l, 17, 3);
	    $resName =~ s/\s//g;
	    my $resNo   = substr($l, 22, 4);
	    $resNo =~ s/\s//g;
	    my $chainId = substr($l, 21, 1);

	    my $aa = "$resName-$resNo-$chainId";
	    if(!exists $uniqAA{$aa}){
#		$fasta .= "\n" if(++$n % 60 == 0 and $fasta ne "");
		$uniqAA{$aa} = $aa;

		if(exists $aa3{$resName}){
		    $fasta .= $aa3{$resName};
		} else {
		    $fasta .= "X";
		}
	    }
	}
    }
    close(F);
    return $fasta;
}
################################################################################
# Needleman-Wunsch algorithm implemented based on Numerical Recipes 3rd Edition, p. 561.  
sub makeNeedlemanWunschAlignment {
    my ($pdbSeq, $uniprotSeq) = @_;
    # scoring scheme
    my $gap      = -4;
    my $skewGap  = 0;

    my $pdbLength = length($pdbSeq);
    my $uniprotLength = length($uniprotSeq);
    
    # initialization
    my @matrix;
    $matrix[0][0] = 0;

    # PDB sequence is in column
    for (my $i = 1; $i <= $pdbLength; $i++) {
	$matrix[$i][0] = $matrix[$i-1][0] + $skewGap;
    }

    # UniProt sequence is in row
    for(my $j = 1; $j <= $uniprotLength; $j++) {
	$matrix[0][$j]  = $matrix[0][$j-1] + $skewGap;
    }

    # determine best scores for each matrix cell
    for (my $i = 1; $i <= $pdbLength; $i++) {
	my $letter1 = substr($pdbSeq, $i-1, 1);
	for(my $j = 1; $j <= $uniprotLength; $j++) {
	    my $letter2 = substr($uniprotSeq, $j-1, 1);
	    my $upScore   = $matrix[$i-1][$j] + (($j == $uniprotLength) ? $skewGap : $gap);
	    my $leftScore = $matrix[$i][$j-1] + (($i == $pdbLength)     ? $skewGap : $gap);
	    my $diagScore = $matrix[$i-1][$j-1] + $blosum62{$letter1}->{$letter2};
	    $matrix[$i][$j] = max(max($upScore, $leftScore), $diagScore);
	}
    }

    # backtrack 
    my $i = $pdbLength;
    my $j = $uniprotLength;
    my $k = 0;
    my $alignUniprotSeq = "";
    my $alignPdbSeq = "";
    while($i > 0 or $j > 0) {
	my $upScore   = -99999999999;
	my $leftScore = -99999999999;
	my $diagScore = -99999999999;
	
	$upScore = $matrix[$i-1][$j] + (($j == $uniprotLength) ? $skewGap : $gap) if($i > 0);
	$leftScore = $matrix[$i][$j-1] + (($i == $pdbLength) ? $skewGap : $gap) if($j > 0);
	if($i > 0 and $j > 0) {
	    my $letter1 = substr($pdbSeq, $i-1, 1);
	    my $letter2 = substr($uniprotSeq, $j-1, 1);
	    $diagScore = $matrix[$i-1][$j-1] + $blosum62{$letter1}->{$letter2};
	}
	if($diagScore >= max($upScore, $leftScore)) {
	    $alignPdbSeq .= substr($pdbSeq, $i-1, 1);
	    $alignUniprotSeq .= substr($uniprotSeq, $j-1, 1);
	    $i--;
	    $j--;
	} elsif ($upScore > $leftScore) {
	    $alignPdbSeq .= substr($pdbSeq, $i-1, 1);
	    $alignUniprotSeq .= "-";
	    $i--;
	} else {
	    $alignUniprotSeq .= substr($uniprotSeq, $j-1, 1);  
	    $alignPdbSeq .= "-";
	    $j--;
	}
    }

    # backtrack creates reverse sequence, so reverse back to normal sequence.
    $alignUniprotSeq = reverse($alignUniprotSeq);
    $alignPdbSeq = reverse($alignPdbSeq);
    print STDERR "\tSequence alignment:\n\t$alignUniprotSeq\n\t$alignPdbSeq\n" if($verbose);
    return &readAlignment($alignPdbSeq, $alignUniprotSeq);
}
################################################################################
sub readAlignment {
    (my $alignUniprotSeq, my $alignPdbSeq) = @_; 

    if(length($alignUniprotSeq) != length($alignPdbSeq)){
	print STDERR "ERROR while reading alignment\n\n";
	exit;
    }

    my $ident = 0;
    my $aln2uniprot;
    my $uniprot2aln;
    my $aln2pdb;
    my $pdb2aln;
    my $aaNo1=1;
    my $aaNo2=1;
    for(my $i = 0; $i < length($alignUniprotSeq); $i++) {
	$ident++ if(substr($alignUniprotSeq, $i, 1) eq substr($alignPdbSeq, $i, 1));

	if(substr($alignUniprotSeq, $i, 1) ne "-"){
	    $aln2uniprot->{$i+1} = $aaNo1;
	    $uniprot2aln->{$aaNo1} = $i+1;
	    $aaNo1++;
	}
	
	if(substr($alignPdbSeq, $i, 1) ne "-"){
	    $aln2pdb->{$i+1} = $aaNo2;
	    $pdb2aln->{$aaNo2} = $i+1;
	    $aaNo2++; 
	}

    }
    $ident /= length($alignUniprotSeq);
    return($ident, $aln2uniprot, $uniprot2aln, $aln2pdb, $pdb2aln);
}
################################################################################
sub pdbResNumbering(){
    (my $pdbFile) = @_; 

    return unless (-e $pdbFile);

    my @a = split(/\n/,`grep ^ATOM $pdbFile`) if(-e $pdbFile);
    my %uniqAA;
    my %map;
    foreach my $l (@a){
	if($l=~/^ATOM/){
	    my $resName = substr($l, 17, 3);
	    $resName =~ s/\s//g;
	    my $resNo   = substr($l, 22, 4);
	    $resNo =~ s/\s//g;
	    my $chainId = substr($l, 21, 1);

	    my $aa = "$resName-$resNo-$chainId";
	    if(!exists $uniqAA{$aa}){
		$uniqAA{$aa} = $aa;
		my $i = keys %uniqAA; 
		$map{$i} = $resNo.$chainId;
	    }
	}
    }
    return \%map;
}
################################################################################
sub pdbSeqNumber2UniprotSeqNumber {
    (my $uniprotSeq, my $pdbSeq) = @_;

    (my $ident, my $aln2uniprot, my $uniprot2aln, my $aln2pdb, my $pdb2aln) =
	&makeNeedlemanWunschAlignment($uniprotSeq, $pdbSeq);
    
    return ($uniprot2aln, $aln2pdb); 
}

##############################################################################
sub fullPath {
    (my $file) = @_;

    my $here = getcwd();
    
    my $dir = `dirname $file | tr -d "\n"`;
    chdir($dir);
    my $fullPath = getcwd()."/".`basename $file | tr -d "\n"`;
    chdir($here);
    return $fullPath;
}
##############################################################################
sub mapPeptide {

    (my $uniprotSeq, my $peptideSeq) = @_;

    my @indices = ();
    if($uniprotSeq =~ /$peptideSeq/i){	
	while ($uniprotSeq =~ /$peptideSeq/gi) {
	    push @indices, [ $-[0], $+[0] ];
	}
    }
    for(my $i=0; $i<=$#indices; $i++){
	$indices[$i][0]++;
    }

    return \@indices;
}

##############################################################################
sub adjustIndices {
    (my $indices, my $trypticity, my $uniprotSeqLength) = @_;

    if($trypticity =~ /N/){
	$indices->[0] = $indices->[1] - 5;
	$indices->[1] = $indices->[1] + 5;
    } elsif($trypticity eq "C"){
	$indices->[1] = $indices->[0] + 5;
	$indices->[0] = $indices->[0] - 5;
    }

    if($indices->[0] < 0){
	$indices->[0] = 0;
    }
    if($indices->[1] > $uniprotSeqLength){
	$indices->[1] = $uniprotSeqLength;
    }
}
##############################################################################
sub adjustPeptide {
    (my $peptides, my $uniprotSeq) = @_;

    my $map;
    foreach my $peptideSeq (sort keys %$peptides){
	my @indices = @{&mapPeptide($uniprotSeq, $peptideSeq)};

	if($#indices > 0) {
	    print STDERR "WARNING: $uniprotSeq has multiple instances of $peptideSeq. ".
		         "Will map only the first instance.\n";
	}

	my $tryptic = "";
	$tryptic .= "N" if(&isNtermTryptic($peptideSeq, substr($uniprotSeq, $indices[0][0]-2,1)));
	$tryptic .= "C" if(&isCtermTryptic($peptideSeq, substr($uniprotSeq, $indices[0][1],1)));

	if($tryptic =~ /N/){
	    $indices[0][0] = $indices[0][1] - 5;
	    $indices[0][1] = $indices[0][1] + 5;
	} elsif($tryptic eq "C"){
	    $indices[0][1] = $indices[0][0] + 5;
	    $indices[0][0] = $indices[0][0] - 5;
	}
	
	if($indices[0][0] < 0){
	    $indices[0][0] = 0;
	}
	if($indices[0][1] > length($uniprotSeq)){
	    $indices[0][1] = length($uniprotSeq);
	}
	my $mapPeptideSeq = substr($uniprotSeq, $indices[0][0], $indices[0][1]-$indices[0][0]+1);
	$map->{$mapPeptideSeq} = "";
    }
    return $map; 
}

##############################################################################
sub pdbResidueNumber {
    (my $uniprotSeqResidueNumber, my $uniprot2aln, my $aln2PDB, my $mapPDBnum) = @_;

    my $alnPos = $uniprot2aln->{$uniprotSeqResidueNumber};

    my $pdbResNo;
    if(defined $alnPos){
	if(exists $aln2PDB->{$alnPos}){
	    $pdbResNo = $mapPDBnum->{$aln2PDB->{$alnPos}};
	}
    }
    return $pdbResNo;
}
##############################################################################
sub pdbTotalResidueNumber {
    (my $uniprotSeqResidueNumber, my $stretchLength, my $uniprot2aln, my $aln2PDB, my $cTerminal) = @_;

    my $alnPos = $uniprot2aln->{$uniprotSeqResidueNumber};

    if(defined $alnPos){
	my $i = $alnPos;
	while($stretchLength-- > 0){
	    if(exists $aln2PDB->{$i}){
		return $aln2PDB->{$i};
	    }
	    if($cTerminal){
		$i++;
	    } else {
		$i--;
	    }
	}
    }
    return;
}
##############################################################################
sub pdbPeptideSeq {

    (my $pdbSeq, my $uniprot2aln, my $aln2PDB, my $start, my $end) = @_;
    my $peptidePDBseq = "";
    for(my $i=$start; $i<=$end; $i++){
	my $pdbTotalResNo = &pdbTotalResidueNumber($i, 1, $uniprot2aln, $aln2PDB, 0);
	if(defined $pdbTotalResNo){
	    $peptidePDBseq .= substr($pdbSeq, $pdbTotalResNo-1, 1);
	} else {
	    $peptidePDBseq .= "-";
	}
    }
    return $peptidePDBseq;
}
##############################################################################
sub printPyMolHeader {
    
    (my $pdbFile, my $uniprotId, my $outDir, my $pdbResNo) = @_;

    return if(!-e $pdbFile);

    if(!-d $outDir){
	print STDERR "Creating directory $outDir\n" if(defined $verbose);
	make_path($outDir);
    }

    my $pmlFile = $pdbFile;
    $pmlFile =~ s/.*\///;
    $pmlFile =~ s/\..*//;
    $pmlFile = "$outDir\/".$uniprotId."_"."$pmlFile.pml";

    my $chainId = "";
    $chainId = substr($pdbResNo, length($pdbResNo)-1, 1) if(defined $pdbResNo);

    open(O, ">$pmlFile") or die "ERROR: Failed to create .pml file $pmlFile: $!.\n";
    print  O "load $pdbFile\n";
    print  O "color grey30\n";
    print  O "hide; set cartoon_fancy_helices, 1; show cartoon; remove solvent; create lig, hetatm; show sticks, lig; show spheres, lig; set sphere_scale, 0.3, lig; cmd.color(\"grey70\", \"lig\"); util.cnc(\"lig\"); zoom;\n";
    printf O "hide everything, all and not chain %s\n", $chainId if($chainId ne "");
    print  O "bg white\n";
    print  O "cmd.alias(\"highlight\", \"";
    close(O);
}
##############################################################################
sub printPyMolScript {
    
    (my $pdbFile, my $uniprotId, my $outDir, my $no, my $uniprotPeptideSeq, my $start, my $end,
     my $uniprot2aln, my $aln2PDB, my $mapPDBnum,  my $tryptic) = @_;

    return if(!-e $pdbFile);

    if(!-d $outDir){
	print STDERR "Creating directory $outDir\n" if(defined $verbose);
	make_path($outDir);
    }

    my $pmlFile = $pdbFile;
    $pmlFile =~ s/.*\///;
    $pmlFile =~ s/\..*//;
    $pmlFile = "$outDir\/".$uniprotId."_"."$pmlFile.pml";

    my $pdbResNo1 = &pdbResidueNumber($start,$uniprot2aln,$aln2PDB,$mapPDBnum);
    my $pdbResNo2 = &pdbResidueNumber($end,$uniprot2aln,$aln2PDB,$mapPDBnum);
    if(!defined $pdbResNo1){
	for(my $i=$start+1; $i<=$end and !defined $pdbResNo1; $i++){
	    $pdbResNo1 = &pdbResidueNumber($i,$uniprot2aln,$aln2PDB,$mapPDBnum);
	}
	if(!defined $pdbResNo1){
	    $pdbResNo1 = "-";
	}
    }
    if(!defined $pdbResNo2){
	for(my $i=$end-1; $i>=$start and !defined $pdbResNo2; $i--){
	    $pdbResNo2 = &pdbResidueNumber($i,$uniprot2aln,$aln2PDB,$mapPDBnum);
	}
	if(!defined $pdbResNo2){
	    $pdbResNo2 = "-";
	}
    }

    my $resNo1 = -1;
    $resNo1 = substr($pdbResNo1, 0, length($pdbResNo1)-1) if(defined $pdbResNo1 and $pdbResNo1 ne "-");
    my $resNo2 = -1;
    $resNo2 = substr($pdbResNo2, 0, length($pdbResNo2)-1) if(defined $pdbResNo2 and $pdbResNo2 ne "-");
    my $chainId = "";
    $chainId = substr($pdbResNo1, length($pdbResNo1)-1, 1) if(defined $pdbResNo1 and $pdbResNo1 ne "-");
    return 0 if($resNo1 == -1 or $resNo2 == -1 or $chainId eq "");

    open(O, ">>$pmlFile") or die "ERROR: Failed to extend .pml file $pmlFile: $!.\n";
    printf O ("select %i-%s-%i-%i, resi %i-%i", $no, $uniprotPeptideSeq, $resNo1, $resNo2, $resNo1, $resNo2);
    printf O (" and chain %s", $chainId) if($chainId ne "");
    printf O ("\\n");
    printf O ("disable %i-%s-%i-%i\\n", $no, $uniprotPeptideSeq, $resNo1, $resNo2);
    printf O ("color yellow, %i-%s-%i-%i\\n", $no, $uniprotPeptideSeq, $resNo1, $resNo2);
    printf O ("label %i-%s-%i-%i and resi %i and chain %s and name CA, %i\\n", $no, $uniprotPeptideSeq, $resNo1, $resNo2, $resNo1, $chainId, $no);
    if($tryptic !~ /N/){
	$nontrypticSelection .= sprintf("select $nontrypticName, (resi %i and chain %s) or nontryptic\ncolor red, $nontrypticName\n", $resNo1, $chainId, $resNo1, $chainId);
    } elsif($tryptic !~ /C/) {
	$nontrypticSelection .= sprintf("select $nontrypticName, (resi %i and chain %s) or nontryptic\ncolor red, $nontrypticName\n", $resNo2, $chainId, $resNo2, $chainId);
    }
    close(O);
    return 1;
}
##############################################################################
sub printPyMolTail {
    (my $pdbFile, my $uniprotId, my $outDir) = @_;

    return if(!-e $pdbFile);

    if(!-d $outDir){
	print STDERR "Creating directory $outDir\n" if(defined $verbose);
	make_path($outDir);
    }

    my $pmlFile = $pdbFile;
    $pmlFile =~ s/.*\///;
    $pmlFile =~ s/\..*//;
    $pmlFile = "$outDir\/".$uniprotId."_"."$pmlFile.pml";

    open(O, ">>$pmlFile") or die "ERROR: Failed to finish .pml file $pmlFile: $!.\n";
    print O "\")\nhighlight\n";
    print O "select $nontrypticName, none\n".$nontrypticSelection."disable $nontrypticName\n" if($nontrypticSelection ne "");
    close(O);
}
###############################################################################
sub isNtermTryptic {
###############################################################################
    # Assesses whether a peptide is tryptic at its Nterm
    # 1st: peptide sequence
    # 2nd: -1 N-terminal AA capital letter
    # return: 1 if N-terminus is tryptic, 0 otherwise.

    (my $peptide, my $nterm) = @_;

    return 1
	if($peptide !~ /^P/
	   and
	   ($nterm =~ /[RK]/
	    or
	    $nterm eq ""));

    return 0;      
}
###############################################################################
sub isCtermTryptic {
###############################################################################
    # Assesses whether a peptide is tryptic at its Cterm
    # 1st: peptide sequence
    # 2nd: +1 C-terminal AA capital letter
    # return: 1 if C-terminus is tryptic, 0 otherwise.

    (my $peptide, my $cterm) = @_;

    return 1
	if($peptide =~ /[RK]$/
	   and
	   ($cterm ne "P"
	    or
	    $cterm eq ""));
    
    return 0;      
}
##############################################################################
##############################################################################
### MAIN
##############################################################################
##############################################################################

&checkInputParameter();
my $peptides = &readInputFile();

foreach my $uniprotId (sort keys %$peptides){
    print STDERR "\n=========================\n".
	           "Processing $uniprotId ...\n" if(defined $verbose);
    
    #############################
    # download UniProt data files
    my $uniprotFile = &downloadUniProt($uniprotId, $fileEndUniProt, $uniprotDir);

    ##############################
    # download UniProt FASTA files
    my $fastaFile = &downloadUniProt($uniprotId, $fileEndFasta, $fastaDir);

    my $cmd = 'grep -v "^>" '.$fastaFile.' | tr -d "\n"';
    my $uniprotSeq = `$cmd`;

    ############
    # download structure files
    print STDERR "Extracting PDB IDs from $uniprotId$fileEndUniProt ...\n" if(defined $verbose);
    my $pdbIds = &pdbIdsFromUniProt("$uniprotDir/$uniprotId$fileEndUniProt");
    print STDERR "FOUND ".(keys %$pdbIds)." PDB chains\n" if(defined $verbose);
    $pdbIds->{999}->{model} = "" if((keys %$pdbIds) == 0);

    foreach my $res (sort keys %$pdbIds){
	foreach my $pdbId (sort keys %{$pdbIds->{$res}}){
	    #####################
	    # download structures
	    last if($maxDownload-- < 0);
	    print STDERR "Downloading structural data for $uniprotId ...\n" if(defined $verbose);
	    my $pdbFile = "";
	    if($pdbId ne "model") {
		foreach my $chainId (sort keys %{$pdbIds->{$res}->{$pdbId}}){
		    print STDERR "Downloading chain ID $chainId of PDB file $pdbId to $pdbDir...\n"
			if(defined $verbose);
		    &downloadPDB($pdbId, $pdbDir, $chainId, 0);
		    $pdbFile = "$pdbDir/$pdbId$chainId$fileEndPdb";
		}
	    }
	    else {
		if(&downloadModBase($uniprotId, $modelDir)){
		    $pdbId = "modbase_".$uniprotId;
		    $pdbFile = "$modelDir/$pdbId$fileEndPdb";
		} else {
		    $pdbId = "";
		    $pdbFile = "";
		}
	    }
	    
	    if($pdbFile eq "") {
		print STDERR "Skipping $uniprotId as neither PDB file nor homology model exists\n\n";
		next;
	    }

	    #########################################
	    # store PDB fasta file in fasta directory
	    my $pdbFileName = $pdbFile;
	    $pdbFileName =~ s/.*\///;
	    $pdbFileName =~ s/\..*//;
	    my $pdbSeq = &pdb2fasta($pdbFile);
	    open (O,">$fastaDir/$pdbFileName.fasta") or die "ERROR: Failed to create fasta file $pdbFileName.fasta: $!.\n";
	    print O ">$pdbId\n$pdbSeq\n";
	    close(O);

	    #########################################################################
	    # determine general residue number mapping between UniProt FASTA sequence and PDB
	    print STDERR "Identifying UniProt <-> PDB residue number mapping ...\n" if(defined $verbose);
	    my $mapPDBnum = &pdbResNumbering($pdbFile);
	    (my $uniprot2aln, my $aln2PDB) = &pdbSeqNumber2UniprotSeqNumber($uniprotSeq, $pdbSeq);

	    ################################################################################
	    # start mapping of peptides on protein structures and extracting all information
	    my $mappedPeptideNo = 0;
	    my $peptideNo = 0;
	    foreach my $uniprotPeptideSeq (sort {$peptides->{$uniprotId}->{$a} <=> $peptides->{$uniprotId}->{$b}} keys %{$peptides->{$uniprotId}}){
		my $peptideSeq = $uniprotPeptideSeq;
		$peptideSeq =~ s/\[.*?\]//g;
		
		$peptideNo++;

		my @indices = @{&mapPeptide($uniprotSeq, $peptideSeq)};
		if($#indices > -1){
		    for(my $i=0; $i<=$#indices; $i++){

			my $tryptic = "";
			$tryptic = "N" if(&isNtermTryptic($peptideSeq, substr($uniprotSeq, $indices[$i][0]-2,1)));
			$tryptic .= "C" if(&isCtermTryptic($peptideSeq, substr($uniprotSeq, $indices[$i][1],1)));
			$tryptic = "-" if($tryptic eq "");

			my $peptideUniprotSeq = substr($uniprotSeq, $indices[$i][0]-1, $indices[$i][1]-$indices[$i][0]+1);
			
			my $pdbPeptideSeq = &pdbPeptideSeq($pdbSeq, $uniprot2aln, $aln2PDB, $indices[$i][0], $indices[$i][1]);
		    
			if($pdbPeptideSeq =~ /[A-Za-z]/) {
			    print STDERR "Mapping peptide $uniprotPeptideSeq onto the protein structure $pdbFile ...\n" if(defined $verbose);
			    &printPyMolHeader($pdbFile, $uniprotId, $pmlDir, &pdbResidueNumber($indices[$i][0], $uniprot2aln, $aln2PDB, $mapPDBnum)) if($mappedPeptideNo == 0);
			    &printPyMolScript($pdbFile, $uniprotId, $pmlDir, $peptideNo, $uniprotPeptideSeq, 
					      $indices[$i][0], $indices[$i][1], $uniprot2aln, $aln2PDB, $mapPDBnum,
					      $tryptic);
			    $mappedPeptideNo++;
			} else {
			    print STDERR "WARNING: None of the aminoacids in the peptide $uniprotPeptideSeq has coordinates in the structure file $pdbFile ...\n" if(defined $verbose);
			}
		    }
		} else {
		    print STDERR "WARNING: $peptideSeq could not be found in $uniprotId. Continuing with next peptide.\n"
			if(defined $verbose);
		}
	    }
	    &printPyMolTail($pdbFile, $uniprotId, $pmlDir) if($mappedPeptideNo > 0);
	}
    }
}
