#!/usr/bin/perl
########################################
##
## Who: Holger Klein & Sergi Sayols
## When: 2017-07-01
## What: custom fastq processor. Count from an input fastq file, how many times
##   a certain pattern occurs.
## 
## @param infile Input fastq file
## @param outdir Output directory
## @param outprefix Prefix output file names will start with
## @param config JSON file with the regular expressions to look for
##
## **NOTE** All arguments are mandatory
##
## @example
##
##   $ ./fastq_processor_riboswitch.pl --infile f1.fastq.gz --outdir ./results --outprefix f1 --config patterns.json
##
## full sequence (Project 531):
## 5’CATTGCAGCGTATTCCCAGTCCCAAACAAACAAAGGCGCGTCCTGGATTCNNNNNAAAACATACCAGATTTCGATCTGGAGAGGTGAAGAATACGACCACCTNNTACATCCAGCTGATGAGTCCCAAATAGGACGAAACGCGCTCAAACAAACAAAGCCTGGTGAAATTGTTATCCGCT -3’
##
## search pattern: 
## - 5nt upstream variable region of 5 + 5nt downstream
## - 5nt downstream of variable region of 2 + 5nt downstream
## - Spacer of (47 - 5 - 5) = 37 nt in between
## GATTC(NNNNN)AAAACATACCAGATTTCGATCTGGAGAGGTGAAGAATACGACCACCT(NN)TACAT
## GATTC(.....)AAAAC.....................................CACCT(..)TACAT
## and reverse complement
## ATGTA(NN)AGGTGGTCGTATTCTTCACCTCTCCAGATCGAAATCTGGTATGTTTT(NNNNN)GAATC
## ATGTA(..)AGGTG.....................................GTTTT(.....)GAATC
##
########################################
use warnings;
use strict;
use Getopt::Long;

binmode STDOUT, ":utf8";
use utf8;

use JSON;

my $infile     = "";
my $outdir     = "";
my $outprefix  = "";
my $configfile = "";
my $help       = 0;

##
## Get comman line arguments
##
sub usage {
    my $message = $_[0];
    my $command = $0;   # get program name from command line
    $command =~ s#^.*/##;
    
    print STDERR (
        $message, 
        "Usage: $command [options]\n" . 
        "Options:\n" . 
        "  --infile f1.fastq.gz\n" . 
        "  --outdir ./results\n" . 
        "  --outprefix f1\n" . 
        "  --config patterns.json\n" . 
        "  --help this help screen"
    );
    
    die("\n")
}

GetOptions('infile=s' => \$infile,
           'outdir=s' => \$outdir,
           'outprefix=s' => \$outprefix,
           'config=s' => \$configfile,
           'help'     => \$help)
           or usage("Invalid commmand line options.\n");

# validate arguments
usage("") if($help);
usage("Incompatible options.\n") if($help && ($infile || $outdir || $outprefix || $configfile));
usage("All arguments are mandatory.\n") if ($infile eq "" || $outdir eq "" || $outprefix eq "" || $configfile eq "");
die("outdir dir $outdir does not exist") unless (-d $outdir);

##
## open input and output files
##
if ($infile =~ /.gz$/) {
    open(IN, "gunzip -c $infile |") || die "Can't open pipe to $infile";
}
else {
    open(IN, $infile) || die "can’t open $infile";
}

open(OUTFILE,     ">" . $outdir . "/" . $outprefix . ".txt");
open(OUTFILEFULL, ">" . $outdir . "/" . $outprefix . "_full.txt");
open(OUTSTATS,    ">" . $outdir . "/" . $outprefix . "_stats.txt");
open(UNMATCHED,   ">" . $outdir . "/" . $outprefix . "_unmatched.txt");

##
## read pattterns from config file
##
my $json;
{
  local $/; #Enable 'slurp' mode
  open my $fh, "<", $configfile or die "File $configfile not found";
  $json = <$fh>;
  close $fh;
}
my $config = decode_json($json);

##
## functions definition
##
sub revcomp {
    my $dna = shift;
    my $rc = reverse($dna);        # reverse the DNA sequence
    $rc =~ tr/ACGTacgt/TGCAtgca/;  # complement the reversed DNA sequence
    return $rc;
}

##
## program body
##
# initialize match counters
my $n_total_match = 0;
my $n_unmatched   = 0;
my %n_fwd_match;
my %n_rvc_match;

foreach (@{$config->{'Patterns'}}) {
    $n_fwd_match{$_->{'Name'}} = 0;
    $n_rvc_match{$_->{'Name'}} = 0;
}

my $matched_pattern;  # name of the pattern matched
my $VRS;              # first switch sequence
my $VRL;              # second switch sequence

# loop through input
while(my $l=<IN>) {

    chomp $l;
    my $match = 0; # is the input pattern in the sequence?
    
    # skip fastq headers and quality values
    next if($l =~ m/^@/);  # header
    next if($l =~ m/^\+/); # spacer between sequence and quality
    next if($l =~ m/[\da-zBDEFHIJKLMOPQRSUVWXYZ]/); # quality string (Phred33)
 
    # loop over the patterns and match the sequence
    foreach (@{$config->{'Patterns'}}) {
        # match on FWD
        if($l =~ m/$_->{'RE'}/) {
            $n_fwd_match{$_->{'Name'}}++;
            $matched_pattern=$_->{'Name'};
            if($_->{'Name'} !~ m/^Spikein/) { # spikein matches go only to global stats
                $match++;
                $VRL=$1;
                $VRS=$2;
            }
            last;  # report only the first matching pattern. Thus, order in the config file matters
        }
        # match on REV
        elsif(revcomp($l) =~ m/$_->{'RE'}/) {
            $n_rvc_match{$_->{'Name'}}++;
            $matched_pattern=$_->{'Name'};
            if($_->{'Name'} !~ m/^Spikein/) { # spikein matches go only to global stats
                $match++;
                $VRL=$1;
                $VRS=$2;
            }
            last;  # report only the first matching pattern. Thus, order in the config file matters
        }
    }
 
    # report the switches sequence
    if($match) {
        print OUTFILE $matched_pattern . "\t" . $VRL . "\t" . $VRS . "\n";                  # without original seq
        print OUTFILEFULL $matched_pattern . "\t" . $VRL . "\t" . $VRS . "\t" . $l . "\n";  # with original seq
        $n_total_match++;
    } 
    else {
        print UNMATCHED $l . "\n";
        $n_unmatched++;
    }
}

# print summary stats
foreach (keys(%n_fwd_match)) {
    print OUTSTATS $_ . " forward:\t" . $n_fwd_match{$_} . "\n";
    print OUTSTATS $_ . " reverse:\t" . $n_rvc_match{$_} . "\n";
}
print OUTSTATS "total:\t". $n_total_match . "\n";
print OUTSTATS "unmatched:\t" . $n_unmatched . "\n";

# and close the file handles
close IN;
close OUTFILE;
close OUTFILEFULL;
close OUTSTATS;
close UNMATCHED;

print "job done!\n";
exit 0;

