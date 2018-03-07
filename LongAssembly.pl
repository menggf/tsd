#!/usr/bin/perl -w
use strict;
#########################################


if(@ARGV==0){
	usage();
	exit;
}

my $output_dir;
my $seq_file;
my $min_fragment_length;
my $insert_seq;
my $genome_ref;
my $cores;
my $min_read;
my $ranking_type;
my $min_merge_overlap=3;
my $merge=0;
my $cutoff_similarity=0.8

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
	elsif ( $ARGV[$i] eq '-r' ) {
		$min_read = $ARGV[ ++$i ];
	}
	elsif ( $ARGV[$i] eq '-o' ) {
		$min_merge_overlap = $ARGV[ ++$i ];
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
$insert_seq="na" if(!defined($insert_seq));

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
if(!-e $output_dir){
  mkdir($output_dir);
}

$genome_ref="/home/gmeng4/genome/hg38/hg38";
system("perl fragment_assembly.pl $output_dir $seq_file $genome_ref $cores $min_fragment_length $insert_seq");
system("perl assembler2.pl $output_dir $ranking_type $min_read $min_fragment_length $min_merge_overlap  $merge $cutoff_similarity");


exit(0);
sub usage{
	print <<'USAGE';
    Command: perl LongAssembly.pl -d output_dir -s seq.fq -l min_fragment_length
                  -i insert_seq.fa -G genome_ref_prefix -p cores -f overlap_score
                  -r min_reads_num -R ranking_type -o min_merge_overlap -m -h
				  
    LongAssembly is a de novo assembly tool for long reads, e.g. PacBio sequencing data.
    It is designed for the sequences with complex structure, e.g. the virus integrated 
    sequences.
	
    Usage:
        -d: the output directory (default: .)
        -s: the sequencing files in fastq format
        -G: the bwa genome reference prefix;
        -l: the minimum fragment length for bwa alignment (default: 200)
        -f: the minimum score for fragment similarity (default: 0.8)
        -i: the extra sequences, e.g. the virus sequences integrated in the genome. (Optional)
        -p: the threads number for parallel computation
        -r: the minimum reads number for final output
        -R: the ranking type, either "length" or "fragment" (default: fragment)
		-o: the minimum overlapped fragment for merging the reads
        -m: merge the reads or not?
        -h: helps
		
    Contact: Guofeng Meng(menggf@gmail.com)
USAGE
}
