####################################################################################
#  alignment.pm @ TSD
#  Info: First-round alignment of PacBio reads.
#  Output: it is the location iformation of pacBio reads after first-round alignment
#
#  Copyright (c) 2017- Guofeng Meng
#
#  EMAIL: menggf@gmail.com



#!/usr/bin/perl -w
use strict;
use TSD::firstround;
package TSD::alignment;

sub new{
	my $class=shift;
	my ($output_dir, $seq_file, $genome_ref, $bwa, $cores, $min_fragment_length, $insert_seq)=@_;
  my $self={};
  srand();
  if(!-e "$output_dir"){
  	mkdir("$output_dir");
  }
  if(!-e "$output_dir/temp_files"){
  	mkdir("$output_dir/temp_files");
  }
  if(!-e "$output_dir/seq_pic"){
  	mkdir("$output_dir/seq_pic");
  }
  if(!-e "$output_dir/insert_ref" and $insert_seq ne "na"){
  	mkdir("$output_dir/insert_ref");
  }
  
  if(!-e "$output_dir/aln_genome.sam"){
    	my $cmd="$bwa mem -t $cores -x pacbio -Y -v 3 -o $output_dir/aln_genome.sam $genome_ref $seq_file";
    	system($cmd)==0 or die;
  }
  
  my $insert_ref="$output_dir/insert_ref/insert";  
  if($insert_seq ne "na"){
    #building the insert reference
    if(!-e "${insert_ref}.bwt"){
    	my $cmd="$bwa index -p $insert_ref $insert_seq";
    	system($cmd)==0 or die;
    }
    if(!-e "$output_dir/aln_insert.sam"){
    	my $cmd="$bwa mem -t $cores -x pacbio -Y -v 3 -o $output_dir/aln_insert.sam $insert_ref $seq_file";
    	system($cmd)==0 or die;
    }
  }
  
  my %range;
  my $fr=TSD::firstround->new($output_dir, $insert_seq, $min_fragment_length);
  my $rg=$fr->{range};
  open OUT,">$output_dir/temp_files/range.txt" or die;
  foreach (keys %$rg){
  	#$range{$_}=$$rg{$_};
  	print OUT "$_\t$$rg{$_}\n";
  }
  undef $rg;
  undef $fr;
  print "Finished first round reading..\n";
  bless($self,$class);
	return($self);
}

sub cigar_len{
	my $s=shift @_;
	my $ll=0;
	my $rr=0;
	my $mm=0;
	my $rmm=0;
	$ll=$1 if($s=~/^(\d+)[HS]/);
	$rr=$1 if($s=~/(\d+)[HS]$/);
	while($s=~/(\d+)(\w)/g){
		if($2 eq "D"  or $2 eq "M" or $2 eq "N" or $2 eq "X" or $2 eq "="){
			$mm += $1;
		}
		if($2 eq "M" or $2 eq "I" or $2 eq "X" or $2 eq "="){
			$rmm += $1;
		}
	}
	return($ll, $mm, $rmm, $rr);
}
sub revs{
	my $s=shift @_;
    $s=reverse($s);
	$s=~tr/ATCGatcg/TAGCtagc/;
	return($s);
}
1;
