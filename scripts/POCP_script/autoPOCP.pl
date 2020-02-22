#!/usr/bin/env perl
use warnings;
use strict;
use Bio::SeqIO;
use File::Spec;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;
use Time::HiRes qw{ sleep };

my $incommand = join( " ", $0, @ARGV );

my $cwd = File::Spec->curdir();
$cwd = File::Spec->rel2abs($cwd);
my $blastout = "$cwd/blast";
my @dbs;
my @chunked;
my ( $size, $pid_cutoff, $coverage );
my $removeNs = 0;    #set to 0 to keep Ns in sequence.
my $help     = 0;
my $man      = 0;
my $finish   = 0;
my $quiet    = 0;
my $log      = 1;
my %idmatch;
my $threads = 1;
my $outfile;
my %matched;
my $blasthitsfile;
my $sge     = 0;
my $blast_v_check = qr/2.2.31|2.[3-9]+.[\d]|3.[\d]+.[\d]+/;
my $outtype = 0;
my $full    = 0;
my $queue   = 'bpp';
my $stime = 0.33;
my $orthofinder = 0;

my ( $svol, $sdir, $sfile ) = File::Spec->splitpath($0);
$sdir .= "scripts/";

my $blastbindir =
  '';    #Set to complete path for the blast bins if not found in PATH
my $program     = $blastbindir . 'blastp';
my $makeblastdb = $blastbindir . 'makeblastdb';
my $scriptdir   = '';
my $blastscript = $scriptdir . 'autoPOCP-blast.pl';
my $elinkpath   = $scriptdir . 'auto_eutil.pl';

if ( !$scriptdir ) {
    $blastscript = $sdir . 'autoPOCP-blast.pl';
    $elinkpath   = $sdir . 'auto_eutil.pl';
}

my $signal = GetOptions(
                         'pid=i'      => \$pid_cutoff,
                         'coverage=i' => \$coverage,
                         'help'       => \$help,
                         'man'        => \$man,
                         'finish'     => \$finish,
                         'quiet'      => \$quiet,
                         'log!'       => \$log,
                         'threads=i'  => \$threads,
                         'outfile=s'  => \$outfile,
                         'sge'        => \$sge,
                         'outtype'    => \$outtype,
                         'full'       => \$full,
                         'queue=s'    => \$queue,
                         'orthofinder' => \$orthofinder
                       );

if ($help) {
    pod2usage( -verbose => 1,
               -output  => ">&STDERR" );
}

if ($man) {
    pod2usage( -verbose => 2,
               -output  => ">&STDERR" );
}

die("Unknown option passed.  Check parameters and try again.\n") if !$signal;

if ( $sge ) {
    my $hostname = `hostname`;
    chomp($hostname);
    if ( $hostname !~ /cgrb\.oregonstate\.local/ ) {
        logger("You do not appear to be on the Oregon State University infrastructure.\n");
        logger("Setting \$sge flag to 0\n");
        $sge = 0;
    } else {
        my @submit_hosts = `qconf -ss`;
        chomp(@submit_hosts);
        my $success = 0;
        foreach my $submit_host (@submit_hosts) {
            if ( $hostname =~ /$submit_host/ ) {
                logger("$hostname is SGE submit host.\n");
                $success = 1;
                last;
            }
        }
        if ( $success == 0 ) {
            logger(
                "$hostname is not submit host.\n"
                . "Retry submission from a submit host below "
                . "(or run the script locally without the -sge flag):\n"
                . join( "\n", @submit_hosts )
                . "\n" 
            );
            exit(-1);
        }
    }   
}

my @infiles;
foreach (@ARGV) {
    my $infile = File::Spec->rel2abs($_);
    push( @infiles, $infile );
}

my %infiles = ();

%infiles = map { $_ => 1 } @infiles if (@infiles);

if ( -s "./queries/queries.txt" ) {
    open my $queries, "<", "./queries/queries.txt"
      or die
      "Unable to open queries.txt (from ./queries/queries.txt). Check permissions and try again.\n";
    my @queries = <$queries>;
    chomp(@queries);
    foreach my $query (@queries) {
        $infiles{$query} = 1;
    }
    close $queries;
}

die("No input specified!") if ( !%infiles );

#Set defaults
$size       //= 1020;    #set for chunk size
$pid_cutoff //= 40;      #set for %ID cutoff
$coverage   //= 50;      #set for coverage cutoff (compared to chunk size)

if ( $pid_cutoff > 70 ) {
    logger(
        "Percent identity cutoff set above 70%. Must be between 0 and 70. Setting to 70.\n"
    );
    $pid_cutoff = 70;
} elsif ( $pid_cutoff < 0 ) {
    logger(
        "Percent identity cutoff set below 0%. Must be between 0 and 70. Setting to 0.\n"
    );
    $pid_cutoff = 0;
}

if ( $coverage > 100 ) {
    logger(
        "Coverage cutoff set to above 100. Must be between 0 and 100. Setting to 100.\n"
    );
} elsif ( $coverage < 0 ) {
    logger(
        "Coverage cutoff set to below 0. Must be between 0 and 100. Setting to 0.\n"
    );
}

my $time = localtime();
logger("Command as submitted at $time:\n$incommand\n");

my @blast_version = `$program -version`;

if ( grep( /$blast_v_check/, @blast_version ) ) {
    logger("Found BLAST version 2.2.31 or greater\n");
} else {
    logger(
        "BLAST version 2.2.31 or higher is REQUIRED!  Make sure the full path is provided or include it in your PATH!\n"
    );
    exit(-1);
}

#Setup blastdbs
my @folders = ( "fasta", "blast", "queries" );

foreach my $folder (@folders) {
    if ( !-d $folder ) {
        my $error = system("mkdir ./$folder");
        die "Unable to generate $folder folder! Check writing permissions. : $!"
          if $error != 0;
    }
}

logger("Saving/Retrieving queries from file\n");

open my $queries, ">", "./queries/queries.txt"
  or die
  "Unable to open queries.txt for writing. Check permissions and try again!\n";
@infiles = sort keys %infiles;
print $queries join( "\n", @infiles ) . "\n";
close $queries;

my %total;

foreach my $infile1 ( sort keys %infiles ) {
    my $proteins = `grep \'>\' $infile1 | wc -l`;
    my ( $vol1, $dir1, $file1 ) = File::Spec->splitpath($infile1);
    $total{$file1} = $proteins;
    foreach my $infile2 ( sort keys %infiles ) {
        next if $infile1 eq $infile2;
        my ( $vol2, $dir2, $file2 ) = File::Spec->splitpath($infile2);
        if ( $file1 eq $file2 ) {
            logger("Found indentically named file in two locations:\n");
            logger("$infile1\n");
            logger("$infile2\n\n");
            logger(
                "Please check these paths and rename one of the files (if appropriate), or remove the duplicate file.\n"
            );
            logger(
                "You must remove the offending file from ./queries/queries.txt to continue.\n"
            );
            exit(-1);
        }
    }
    if ( !-s $infile1 ) {
        logger("Unable to find input file $infile1.\n");
        logger(
            "Check the path name and correct it in ./queries/queries.txt before continuing.\n"
        );
        exit(-1);
    }
}

$time = localtime();

logger("Setting up blast databases at $time\n");

foreach my $infile (@infiles) {
    my ( $volume, $dir, $file ) = File::Spec->splitpath($infile);
    if (    ( -e "db/$file.phr" )
         && ( -e "db/$file.pin" )
         && ( -e "db/$file.psq" ) )
    {
        logger("Blast database already present for db/$file. Skipping...\n");
        push( @dbs, "$cwd/db/$file" );
        next;
    } else {
        my $output = `$makeblastdb -in $infile -dbtype prot -out db/$file`;
        logger( $output . "\n" );
        push( @dbs, "$cwd/db/$file" );
        if ( $? != 0 ) {
            logger("Problem generating blast databases!\n");
            exit(-1);
        }
    }
}

$time = localtime();

logger("Done at $time.\n");

if ( $finish == 0 ) {

    my $outfmt = q{6 std qcovhsp};

    my $pm = Parallel::ForkManager->new( $threads - 1 );

    $time = localtime();
    logger("BLAST searches started at $time\n");

    foreach my $querypath (@infiles) {
        my ( $vol, $dir, $query ) = File::Spec->splitpath($querypath);
        foreach my $dbpath (@dbs) {
            $pm->start and next;
            my ( $vol1, $dir1, $subject ) = File::Spec->splitpath($dbpath);
            if ( $query ne $subject ) {
                my $output = "./blast/";
                $output .= "Blast" . $query . "_vs_" . $subject . ".txt";
                if ( -s $output ) {
                    logger("BLAST output $output already found!\n");
                } else {
                    if ( -e $output ) {
                        `rm -f $output`;
                    }
#                    my @command =
#                      join( " ",
#                            $program,                 "-db",
#                            qq{\\'\'$dbpath\'\\'},    "-query",
#                            qq{\\'\'$querypath\'\\'}, "-outfmt",
#                            $outfmt,                  "-evalue",
#                            "1e-5",                   "-num_threads",
#                            1,                        "-max_target_seqs",
#                            "1",                      "-max_hsps",
#                            "1" );
#                    logger("Running BLAST command:\n@command\n");
#                    @command = join( " ",
#                                     $blastscript, $coverage, $pid_cutoff,
#                                     $output,      $query,    $subject,
#                                     @command );
                    if ( $sge == 0 ) {
                        my @command =       ($program,                 
                                            "-db",              qq{$dbpath},    
                                            "-query",           qq{$querypath}, 
                                            "-outfmt",          qq{$outfmt},                  
                                            "-evalue",          "1e-5",   
                                            "-num_threads",     1,         
                                            "-max_target_seqs", "1",      
                                            "-max_hsps",        "1",      
                                            "-out",             $output);
                        logger("Running BLAST command:\n@command\n");
                        my $check = system(@command);
                        if ( $check != 0 || ! -s $output ) {
                            `rm -f $output`;
                            logger("BLAST search between $query and $subject failed.\n");
                            logger("Check your input files and try again.\n");
                            exit(-1)
                        }

                    } else {
                        #Only works on Oregon State University CGRB Infrastructure
                        my @command =       ($program,                 
                                            "-db",              qq{$dbpath},    
                                            "-query",           qq{$querypath}, 
                                            "-outfmt",          qq{\'$outfmt\'},                  
                                            "-evalue",          "1e-5",   
                                            "-num_threads",     1,         
                                            "-max_target_seqs", "1",      
                                            "-max_hsps",        "1",      
                                            "-out",             $output);
                        logger("Running BLAST command through SGE:\n@command\n");
                        @command = join( " ",
                                         "SGE_Batch",                "-r",
                                         "sge.${query}_vs_$subject", "-q",
                                         $queue,                      "-c",
                                         qq{"@command"},             "-Q" );
                        my @runoutput = `@command`;
                        sleep($stime);
                    }

                }

            }
            $pm->finish;
        }
    }
    $pm->wait_all_children;

    $time = localtime();
    logger("BLAST searches completed at $time\n");

    logger("To complete analysis, use the -finish flag.\n");
    exit(0);
}

#my %results;
#my @output;

my %results;
my %tmpresults;
my %genomes;
my $tmpfile = 'pocp.out.tmp';
my %seenfh;
my %filenames;

if ( -s $tmpfile ) {
    logger("Temp file ($tmpfile) found. Extracting previous results now.\n");
    open my $tmpfh, "<", "$tmpfile" or die "Unable to open $tmpfile : $!\n";
    while(<$tmpfh>) {
        my $line = $_;
        chomp($line);
        my ($query, $subject, $pocp, $hits, $file) = split("\t",$line);
        if ($file) {
            if (-s "./blast/$file") {
                $seenfh{$file} = 1;
            }
        }
        $tmpresults{$query}{$subject}{'pocp'} = $pocp;
        $tmpresults{$query}{$subject}{'hits'} = $hits;
        $genomes{$query} = 1;
    }
}

my @blastout = `find ./blast/ -type f -exec basename {} \\;`;

chomp(@blastout);

my $blastfiles = @blastout;

if ($blastfiles < 2) {
    logger("Unable to find BLAST output in ./blast/\n");
    logger("Check input files and ./queries/queries.txt to ensure proper files are being examined.\n");
    exit(-1);
}

$time = localtime();


logger("Found $blastfiles BLAST output files to parse. Starting at $time.\n");

foreach my $boutfile (sort @blastout) {
    next if $seenfh{$boutfile};

    my @data = blast_summary($boutfile);
    if ($data[0] eq '0') {
        next;
    }
    $results{$data[0]}{$data[1]} = $data[2];
    $genomes{$data[0]} = 1;

}
#open my $tmpfile, "<", "pocp.tmp"
#  or die "Unable to find temp file pocp.tmp. Cannot calculate results!\n";
#
#while (<$tmpfile>) {
#    my $line = $_;
#    chomp($line);
#    my ( $query, $subject, $hits ) = split( "\t", $line );
#    if ( !$query || !$subject || $hits =~ /[A-z]/ || !$hits ) {
#        next;
#    }
#    push( @output, $line );
#    if ( !defined( $total{$query} ) || !defined( $total{$subject} ) ) {
#        next;
#    }
#
#    $results{$query}{$subject} = $hits;
#}
#
#close $tmpfile;

logger("Printing results\n");

my $outfh;
my $hitsfh;

open my $tmpfh, ">>", "$tmpfile" or die "Unable to open $tmpfile : $!\n";

if ($outfile) {
    open $outfh, ">", $outfile
      or die "Unable to open output file $outfile : $!";
} else {
    $outfh = \*STDOUT;
}

my @genomes = sort keys %total;

my @print;

foreach my $genome (@genomes) {
    push( @print, trim($genome) );
}

if ( $outtype == 0 ) {
    print $outfh join( "\t", "", @print ) . "\n";
} else {
    print $outfh join( "\t", "genome1", "genome2", "pocp" ) . "\n";
}


for my $g1 (@genomes) {
    if ( $outtype == 0 ) {
        print $outfh trim($g1) . "\t";
    }
    for ( my $j = 0 ; $j < scalar(@genomes) ; $j++ ) {
        my $g2 = $genomes[$j];
        if ( $g1 eq $g2 ) {
            print $outfh '----------'."\t";
            next;
        } else {
            my $pocp = "NA";
            my $trg1 = (defined($tmpresults{$g1}{$g2}{'pocp'}) ? $tmpresults{$g1}{$g2}{'pocp'} : 0);
            my $trg2 = (defined($tmpresults{$g2}{$g1}{'pocp'}) ? $tmpresults{$g2}{$g1}{'pocp'} : 0);
            if ($trg1 == $trg2 && $trg1 > 0) {
                $pocp = $trg1;
            } else {
                my $print = 0;
                my ($hitsg1, $hitsg2);
                $hitsg1 = (defined($tmpresults{$g1}{$g2}{'hits'}) ? $tmpresults{$g1}{$g2}{'hits'} : $results{$g1}{$g2});
                $hitsg2 = (defined($tmpresults{$g2}{$g1}{'hits'}) ? $tmpresults{$g2}{$g1}{'hits'} : $results{$g2}{$g1});

                $pocp = sprintf( "%.2f", ( $hitsg1 + $hitsg2 ) / 
                          ( $total{$g1} + $total{$g2} ) * 100);
                if (defined($results{$g1}{$g2})) {
                    print $tmpfh join("\t",$g1,$g2,$pocp,$results{$g1}{$g2},$g1."_vs_".$g2.".txt")."\n";
                }
            }

            if ( $outtype == 0 ) {
                print $outfh $pocp . "\t";
            } else {
                print $outfh join( "\t", $g1, $2, $pocp ) . "\n";
            }
        }
    }

    if ( $outtype == 0 ) {
        print $outfh "\n";
    }
}

close $tmpfh;
close $outfh;

$time = localtime();

logger("Finished at $time.\n");

sub logger {
    my $message = shift;
    print STDERR $message unless $quiet == 1;
    if ( $log == 1 ) {
        open my $fh, ">>", "pocp.log" or die "pocp.log is unavailable : $!";
        print $fh $message;
        close $fh;
    }
}

sub trim {
    my $line = shift;

    $line =~ s/\.[^\.]+$//;

    return $line;
}

sub split_names {
    my $name = shift;
    my @fs = split('_vs_',$name);
    my $qu = $fs[0];
    my $su = $fs[1];
    $qu =~ s/^Blast//;
    $su =~ s/\.txt$//;
    return ($qu, $su);
}

sub blast_summary {
    my $blast   = shift;
    my ($query, $subject) = split_names($blast);
    $blast = "./blast/$blast";
    
    if ( ! -e "$blast" ) {
        logger("Unable to find blast file $blast.\n");
        return 0;
    }

    open my $infh, "<", $blast or die "Unable to open input file $blast : $!";
    my $hits = 0;
    while (<$infh>) {
        my $line = $_;
        chomp($line);
        next if $line =~ /#/;
        #qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp

        my @data      = split( "\t", $line );
        my $percentid = $data[2];
        my $qcov      = $data[12];

        if ( $percentid >= $pid_cutoff && $qcov >= $coverage ) {
            $hits++;
        }

    }
    close $infh;
    return ($query, $subject, $hits);
}

__END__

=head1 NAME

pocp.pl - Perform pairwise POCP (Percentage Of Conserved Proteins) comparisons between any number of sequenced genomes

=head1 SYNOPSIS

pocp.pl input[n].fasta input[n-1].fasta ... input[2].fasta input[1].fasta 

=head1 OPTIONS

Defaults shown in square brackets.  Possible values shown in parentheses.

=over 8

=item B<-help|h>

Print a brief help message and exits.

=item B<-man>

Print a verbose help message and exits.

=item B<-quiet>

Turns off progress messages. Use noquiet to turn messages back on.

=item B<-log|nolog> [logging on]

Using -nolog turns off logging. Logfile is in the format of pocp.log.

=item B<-threads> [1]

*Recommended* - Using multiple threads will significantly speed up the BLAST searches.

=item B<-keys>

Re-generate keyfile from NCBI/local FASTA files.  Generally unnecessary unless pocp.keys has been altered/removed.

=item B<-outfile> 

Optional.  Path to output file. By default, output is directed to pocp.out.

=item B<-coverage> [50] (0-100)

Percentage of query coverage cutoff.

=item B<-pid> [40] (0-70)

Percent identity cutoff for BLAST search results.

=back

=cut
