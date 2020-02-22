#!/usr/bin/perl
#This script is strictly to run BLAST searches for Auto ANI
#Generating the data in this manner allows for felxibility
#in parallelization

use warnings;
use strict;

my $coverage = shift;
my $pid_cutoff = shift;
my $output = shift;
my $query = shift;
my $subject = shift;
my @command = join( " ", @ARGV );
my %results;
my %genomes;
my $hits = 0;

@command = join(" ",@command,'|','tee',$output);

my @output = `@command`;

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

open my $tmpout, ">>", "pocp.tmp" or die "Unable to open pocp.tmp for writing!\n";

print $tmpout join("\t",$query,$subject,$hits)."\n";

close $tmpout;

