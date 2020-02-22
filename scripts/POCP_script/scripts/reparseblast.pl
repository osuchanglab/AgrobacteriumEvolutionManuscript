#!/usr/bin/perl
#This script is strictly to run BLAST searches for Auto ANI
#Generating the data in this manner allows for felxibility
#in parallelization

use warnings;
use strict;

my $infile = shift
my $coverage = 50;
my $pid_cutoff = 40;
my %results;
my %genomes;
my $hits = 0;



open my $handle, '<', $infile;
chomp(my @output = <$handle>);
close $handle;


foreach my $line (@output) {
    chomp($line);
    next if ( $line =~ /#/ );
    my @data      = split( "\t", $line );
    my $percentid = $data[2];
    my $qcov      = $data[3];

    if ( $percentid >= $pid_cutoff && $qcov >= $coverage ) {
        $hits++;   
    };

}

open my $tmpout, ">>", "pocpnew.tmp" or die "Unable to open pocpnew.tmp for writing!\n";

print $tmpout join("\t",$query,$subject,$hits)."\n";

close $tmpout;

