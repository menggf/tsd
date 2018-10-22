####################################################################################
#  LongAssembly.pl @ TSD
#  This is the interface file for TSD package.
#
#  Copyright (c) 2017- Guofeng Meng
#
#  EMAIL: menggf@gmail.com


#!/usr/bin/perl -w
use strict;
my $lib; 
BEGIN {
  $lib = "lib"; # the file library path
}
use lib $lib; 
use TSD::alignment;
use TSD::assemblyReads;
use TSD::assemblyCF;

if(@ARGV==0){
	usage();
	exit;
}

my $output_dir;
my $seq_file;
my $min_fragment_length;
my $insert_seq="na";
my $genome_ref;
my $cores=1;
my $min_read;
my $ranking_type;
my $min_merge_overlap=3;
my $merge=0;
my $cutoff_similarity=0.7;
my $bwa=`which bwa`;
chomp $bwa;
if($bwa=~/no bwa/){
  print "Error: Please install bwa or add bwa to system \$PATH\n";
	usage();
	exit;
}

for ( my $i = 0 ; $i <= $#ARGV ; $i++ ) {
	if ( $ARGV[$i] eq '-d' ) {
		$output_dir = $ARGV[ ++$i ];
	}
	elsif ( $ARGV[$i] eq '-s' ) {
		$seq_file = $ARGV[ ++$i ];
	}
	elsif ( $ARGV[$i] eq '-l' ) {
		$min_fragment_length = $ARGV[ ++$i ];
	}
	elsif ( $ARGV[$i] eq '-i' ) {
		$insert_seq = $ARGV[ ++$i ];
	}
	elsif ( $ARGV[$i] eq '-p' ) {
		$cores = $ARGV[ ++$i ];
	}
	elsif ( $ARGV[$i] eq '-G' ) {
		$genome_ref = $ARGV[ ++$i ];
	}
	elsif ( $ARGV[$i] eq '-R' ) {
		$ranking_type = $ARGV[ ++$i ];
	}
	elsif ( $ARGV[$i] eq '-r' ) {
		$min_read = $ARGV[ ++$i ];
	}
	elsif ( $ARGV[$i] eq '-o' ) {
		$min_merge_overlap = $ARGV[ ++$i ];
	}
	elsif ( $ARGV[$i] eq '-e' ) {
		$bwa = $ARGV[ ++$i ];
	}
	elsif ( $ARGV[$i] eq '-f' ) {
		$cutoff_similarity = $ARGV[ ++$i ];
	}
	elsif ( $ARGV[$i] eq '-m' ) {
		$merge = 1;
	}
	elsif ( $ARGV[$i] eq '-h' ) {
		usage();
		exit();
	}
}

$ranking_type = "fragment" if(!defined($ranking_type));
$min_read = 5 if(!defined($min_read));
$min_fragment_length=200 if(!defined($min_fragment_length));
$cores=1 if(!defined($cores));


if (!defined($output_dir)) {
	print "Error: Please specifiy the output directory\n";
	usage();
	exit;
}
if (!defined($seq_file)) {
	print "Error: Please specifiy the fastq sequencing file\n";
	usage();
	exit;
}
if (!defined($genome_ref)) {
	print "Error: Please specifiy the bwa genome reference prefix\n";
	usage();
	exit;
}
if (!-e "$genome_ref.bwt") {
	print "Error: Cannot find bwa genome reference: $genome_ref\n";
	usage();
	exit;
}
if(!-e $output_dir){
  mkdir($output_dir);
}

if(!defined($insert_seq) or  $insert_seq=="na"){
  $insert_seq="na";
}
else{
  open INPUT,$insert_seq or die "Did not find $insert_seq";
  open OUT,">$output_dir/insert.fa" or die;
  while(my $str=<INPUT>){
    if($str=~/^>/){
      print OUT ">insert\n";
    }
    else{
      print OUT $str;
    }
  }
  close INPUT;
  close OUT;
  $insert_seq="$output_dir/insert.fa";
}

#$genome_ref="/home/meng/work/genome/hg38.fa";
TSD::alignment ->new($output_dir, $seq_file, $genome_ref, $bwa, $cores, $min_fragment_length, $insert_seq);
TSD::assemblyReads->new($output_dir, $seq_file, $genome_ref, $bwa, $cores, $min_fragment_length, $insert_seq);
TSD::assemblyCF->new($output_dir, $ranking_type, $min_read, $min_fragment_length, $min_merge_overlap,  $merge, $cutoff_similarity, $lib,$insert_seq);


exit(0);
sub usage{
	print <<'USAGE';
    Command: perl LongAssembly.pl -d output_dir -s seq.fq -l min_fragment_length
                  -i insert_seq.fa -G genome_ref_prefix -p cores -f overlap_score
                  -r min_reads_num -R ranking_type -o min_merge_overlap -e path_bwa
                  -m -h
				  
    LongAssembly is a de novo assembly tool for long reads, e.g. PacBio sequencing data.
    It is designed for the targeted sequences with complex structure, e.g. the virus 
    integrated sequences.
	
    Usage:
        -d: the output directory (default: .)
        -s: the sequencing files in fastq format
        -G: the bwa genome reference prefix;
        -e: the path of bwa
        -l: the minimum fragment length for bwa alignment (default: 200)
        -f: the minimum score for fragment similarity (default: 0.7)
        -i: the extra sequences, e.g. the virus sequences integrated in the genome. (Optional)
        -p: the threads number for parallel computation
        -r: the minimum reads number for final output
        -R: the ranking type, either "length" or "fragment" (default: fragment)
        -o: the minimum overlapped fragment for merging the reads
        -m: merge the reads or not?
        -h: helps
		
    Note: TSD support continuing analysis by automaticly detecting the output files in output 
    directory. For examples, if "aln_genome.sam" exists in output directory, TSD will not do 
    first-round genomic alignment again. Therefore, if users want to do the whole analysis from
    beginning, please delete the output directory completely.
    
    Contact: Guofeng Meng(menggf@gmail.com)
USAGE
}
