#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 09.09.2016

###############################################################################
###############################################################################
### Retrieves the starting and end position of peptides within a list of    ###
### FASTA sequences. In addition, it outputs information whether the        ###
### peptide is tryptic, halftryptic, or non-tryptic.                        ###
###############################################################################
###############################################################################

use strict;
use warnings;
use Getopt::Long;

my (
    # variable for parameters which are read in from commandline
    $help,
    $peptideFile,
    $col,
    $fastaFile,
    $verbose,
    $hasHeader,
   );

##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "help!"     => \$help,        # Print this help
    "in=s"      => \$peptideFile, # Path to a tab delimited text file, which holds in one column a list of peptide sequences
    "col:i"     => \$col,         # Column number of peptide sequences. Numbering starts with 1. Default is 2
    "header!"   => \$hasHeader,   # Ignore header in peptide sequence file
    "fasta=s"   => \$fastaFile,   # Path to a FASTA formatted proteome sequence text file
    "verbose!"  => \$verbose,     # Outputs additional information such as a table header
) or die "\nTry \"$0 -h\" for a complete list of options.\n\n";

##############################################################################
$col = 2 if(!defined $col);
$col--;
##############################################################################
# help
if ($help) {printHelp(); exit}

# Catch user input errors
# Report if required options are missing.
if(!defined $peptideFile or !defined $fastaFile){
    print STDERR "\n";
    print STDERR "Please provide the path to the peptide file via the -in option.\n"
	if(!defined $peptideFile);
    print STDERR "Please provide the path to the FASTA file via the -fasta option.\n"
	if(!defined $fastaFile);
    print STDERR "For more details on the list of options, use the -h option.\n\n";
    exit;
} else {
    # Report if file is missing
    if(! -e $peptideFile or ! -e $fastaFile){
	print STDERR "\n";
	print STDERR "No peptide file found at path $peptideFile\n" unless(-e $peptideFile);
	print STDERR "No FASTA file found at path $fastaFile\n" unless(-e $fastaFile);
	print STDERR "\n";
	exit;
    }
}

##############################################################################
### SUBROUTINES ##############################################################
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

    $usage = "$0 -in myoglobin_lip_peptides.tsv -fasta horse_proteome.fasta\n";

    print "\nRetrieves the starting and end position of peptides within a list of FASTA sequences.\n";
    print "\nOutput format: PepSeq\tProtName\tStartPos\tEndPos\t".
	  "FreqInProt\tFreqAmongProts\tTryptic\tExtendPepSeq\n";
    print "\nUsage:\n" .  $usage . "\n";

    print "Valid options are:\n\n";
    open(MYSELF, "$0") or
      die "Cannot read source code file $0: $!\n";
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

###############################################################################
sub isTryptic {
###############################################################################
    # Assesses whether a peptide is fully, half or non-tryptic
    # param:
    # 1st: peptide sequence
    # 2nd: -1 N-terminal AA capital letter
    # 3rd: +1 C-terminal AA capital letter
    # return: 0 if non-tryptic, 0.5 if half-tryptic, 1 if tryptic

    (my $peptide, my $nterm, my $cterm) = @_;


    if(&isNtermTryptic($peptide, $nterm)
       and
       &isCtermTryptic($peptide, $cterm)){
	return 1;
    }
    if(!&isNtermTryptic($peptide, $nterm)
       and
       !&isCtermTryptic($peptide, $cterm)){
	return 0;
    }
    return 0.5;
}


##############################################################################
### END OF SUBROUTINES########################################################
############################################################################## 

############
### MAIN ###
############

if($fastaFile =~ /\.gz$/) {
    open(F1, "zcat < $fastaFile |") or die "Failed to open fasta file \"$fastaFile\"\n";
} else {
    open(F1, $fastaFile) or die "Failed to open fasta file \"$fastaFile\"\n";
}

# Collect all FASTA sequences from the FASTA file into a hash
my $header = "";
my $seq = "";
my %seqs = ();
while(my $line = <F1>){
    chomp($line);
    if($line =~ /^>/){
	$header = $line;
    } else {
	$seqs{$header} .= $line;
    }	
}
close(F1);    


open(F2, $peptideFile) or die "Failed to open peptide file \"$peptideFile\"\n";
$header = <F2> if(defined $hasHeader);
my %peptideSeqs = ();
while(my $line = <F2>){
    chomp($line);

    my @a = split(/\t/, $line);
    @a = split(/\,/, $line) if($peptideFile =~ /.csv$/);

    my $peptide = $a[$col];
    # in case input file is comma seperated and fields have quotes
    $peptide =~ s/\"//g;
    $peptideSeqs{$peptide} = 1;
}
close(F2);

print "#PepSeq\tProtName\tStartPos\tEndPos\t".
      "FreqInProt\tFreqAmongProts\tTryptic\tExtendPepSeq\n" if(defined $verbose);

# Search for peptides in the FASTA sequence hash
foreach my $peptide (sort keys %peptideSeqs) {
    my @startPos = ();
    my @endPos = ();
    my @proteinNames = ();
    my @tryptic = ();
    my $frequencyAmongProteins;
    my %frequencyWithinProtein = ();
    my @extendedPeptideSeq = ();

    my %n = (); # hash to store the protein names that contain the peptide sequence
    # search each peptide in all protein sequences
    foreach $header (sort keys %seqs){
	my $name = $header;
	$name =~ s/^>//;
	$name =~ s/ .*//; 

	$frequencyWithinProtein{$name} = 0;
	# search peptide within a single protein and store the number of times
	# it occurred within that single protein together with its starting
	# and end position
	my $protein = $seqs{$header};
	while ($protein =~ /(.{0,1})$peptide(.{0,1})/g) {
	    my $nterm = $1;
	    my $cterm = $2;
	    $frequencyWithinProtein{$name}++;

	    push(@startPos, $-[0]);
	    push(@endPos, $+[0]);
	    push(@proteinNames, $name);
	    $n{$name} = "";
	    
	    push(@tryptic, &isTryptic($peptide, $nterm, $cterm)); 

	    # Extract the N and C-terminal AA of the peptide within its
	    # protein sequence, with the dash symbol (-) representing the
	    # N or C terminus respectively. In any case, the N and C-terminal
	    # AA is seperated from the peptide sequence with a dot.
	    $nterm = "-" if($nterm eq "");
	    $cterm = "-" if($cterm eq "");

	    push(@extendedPeptideSeq, "$nterm.$peptide.$cterm");
	}
    }
    # store the number of times the peptide was observed with all proteins
    $frequencyAmongProteins = keys %n;

    # create output
    if(@startPos == 0) {
	print "#$peptide\t-\t-\t-\t-\t-\t-\t-\n";
    } else {
	for(my $i = 0; $i < @startPos; $i++){
	    print "$peptide\t$proteinNames[$i]\t".($startPos[$i]+1)."\t".$endPos[$i]."\t$frequencyWithinProtein{$proteinNames[$i]}\t$frequencyAmongProteins\t$tryptic[$i]\t$extendedPeptideSeq[$i]\n";
	} 
    }
}
close(F2);    

exit;
