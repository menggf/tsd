#!/usr/bin/perl -w
use strict;

package firstround;

sub new{
	my $class=shift;
	my ($output_dir,$insert_seq, $min_fragment_length)=@_;
	my $self={};
	my %all1;
	my %all2;
	my %range;
	my %ids;
	open GENOME, "$output_dir/aln_genome.sam" or die;
	while(my $str=<GENOME>){
		next if($str=~/^@/);
		my @str=split(/\t/,$str);
		next if($str[2] eq "*"  or $str[1] > 128);
		$all1{$str[0]}="$str[0]\t$str[1]\t$str[2]\t$str[3]\t$str[5]\t$str[9]";
		$ids{$str[0]}=1;
		last;
	}
	while(my $str=<GENOME>){
		my @str=split(/\t/,$str);
		next if($str[2] eq "*"  or $str[1] > 128);
		$all1{$str[0]}="$str[0]\t$str[1]\t$str[2]\t$str[3]\t$str[5]\t$str[9]";
		$ids{$str[0]}=1;
	}
	close GENOME;

	if($insert_seq ne "na"){
		open INSERT, "$output_dir/aln_insert.sam" or die;
		while(my $str=<INSERT>){
			next if($str=~/^@/);
			my @str=split(/\t/,$str);
			next if($str[2] eq "*" or $str[1] > 128);
			$all2{$str[0]}="$str[0]\t$str[1]\t$str[2]\t$str[3]\t$str[5]\t$str[9]";
			$ids{$str[0]}=1;
			last;
		}
		while(my $str=<INSERT>){
			my @str=split(/\t/,$str);
			next if($str[2] eq "*" or $str[1] > 128);
			$all2{$str[0]}="$str[0]\t$str[1]\t$str[2]\t$str[3]\t$str[5]\t$str[9]";
			$ids{$str[0]}=1;
		}
		close INSERT;
	}
	print "Finishing reading alignment\n";

	my @ids=keys(%ids);
	undef %ids;

	# 0id / 1tag / 2chr / 3pos / 4cigar / 5seq
	open(TEMP,">$output_dir/temp_files/temp_seq_0.fa") or die;
	foreach my $id (@ids){
		if(defined($all1{$id}) and !defined($all2{$id})){
			#print "1\t$id\n";
			my @str1=split(/\t/,$all1{$id});
			my $strand1=1;
			#$strand1=-1 if($str1[1]==16);
			my ($lefts, $match_len,$match_rlen,$rights)=&cigar_len($str1[4]);
			my $from=$str1[3];
			my $to=$str1[3] + $match_len;
			$range{$id}="$str1[2]\t$strand1\t$from\t$to"; 
			if($lefts > $min_fragment_length){
				my $seq_left=substr($str1[5], 0, $lefts);
				print TEMP ">$str1[0]#left\n$seq_left\n";
			}
			if($rights > $min_fragment_length){
				my $seq_right=substr($str1[5], $lefts + $match_rlen-1 ,$rights);
				print TEMP ">$str1[0]#right\n$seq_right\n";
			}    
		}
		elsif(!defined($all1{$id}) and defined($all2{$id})){
			#print "2\t$id\n";
			my @str2=split(/\t/,$all2{$id});
			my $strand2=1;
			#$strand2=-1 if($str2[1]==16);
			my ($lefts, $match_len,$match_rlen,$rights)=&cigar_len($str2[4]);
			#print "$lefts\t$match_len\t$match_rlen\t$rights\t$str2[4]\t",length($str2[5]),"\n";
			my $from=$str2[3];
			my $to=$str2[3] + $match_len;
			$range{$id}="$str2[2]\t$strand2\t$from\t$to"; 
			if($lefts > $min_fragment_length){
				my $seq_left=substr($str2[5], 0, $lefts);
				print TEMP ">$str2[0]#left\n$seq_left\n";
			}
			if($rights > $min_fragment_length){
				my $seq_right=substr($str2[5], $lefts + $match_rlen-1 ,$rights);
				print TEMP ">$str2[0]#right\n$seq_right\n";
			}   
		}
		elsif(defined($all1{$id}) and defined($all2{$id})){
			#print "3\t$id\n";
			my @str1=split(/\t/,$all1{$id});
			my @str2=split(/\t/,$all2{$id});
			my $strand2=1;
			$strand2=-1 if($str2[1]==16);
			my $strand1=1;
			$strand1=-1 if($str1[1]==16);
			if($strand1==-1 and $strand2==-1){
				$strand1=1;
				$strand2=1;
			}
			elsif($strand1==-1 and $strand2==1){
				$strand1=1;
				$strand2=-1;
			}
			my ($lefts1, $match_len1,$match_rlen1,$rights1)=cigar_len($str1[4]);
			my ($lefts2, $match_len2,$match_rlen2,$rights2)=cigar_len($str2[4]);
			my $from2=$str2[3];
			my $to2=$str2[3] + $match_len2;
			my $from1=$str1[3];
			my $to1=$str1[3] + $match_len1;	
			if($strand2 ==-1){
				my $tmp=$lefts2;
				$lefts2=$rights2;
				$rights2=$tmp;
			}
			if($lefts1 <= $lefts2){
				if($lefts1 > $min_fragment_length){
					my $seq_right=substr($str1[5], 0 ,$lefts1);
					print TEMP ">$str1[0]#left\n$seq_right\n";
				}
				my $gap=$lefts2 -$lefts1 - $match_rlen1;
				$range{$str2[0]}="$str1[2]\t$strand1\t$from1\t$to1\t1000\t$gap\t$str2[2]\t$strand2\t$from2\t$to2"; 
				if( $gap > $min_fragment_length){
					my $seq_mid=substr($str1[5], $lefts1 + $match_rlen1-1, $gap);
					print TEMP ">$str2[0]#1000\n$seq_mid\n";
				}
				if($rights2 > $min_fragment_length){
					my $seq_right=substr($str1[5], $lefts2 +$match_rlen2-1 ,$rights2);
					print TEMP ">$str1[0]#right\n$seq_right\n";
				}
			}
			if($lefts1 > $lefts2){
				if($lefts2 > $min_fragment_length){
					my $seq_right=substr($str1[5], 0 ,$lefts2);
					print TEMP ">$str1[0]#left\n$seq_right\n";
				}
				my $gap=$lefts1 -$lefts2- $match_rlen2;
				$range{$str2[0]}="$str2[2]\t$strand2\t$from2\t$to2\t1000\t$gap\t$str1[2]\t$strand1\t$from1\t$to1"; 
				if( $gap > $min_fragment_length){
					my $seq_mid=substr($str1[5], $lefts2 + $match_rlen2-1 , $gap);
					print TEMP ">$str2[0]#1000\n$seq_mid\n";
				}
				if( $rights1 > $min_fragment_length){
					my $seq_right=substr($str1[5], $lefts1 +$match_rlen1-1 ,$rights1);
					print TEMP ">$str1[0]#right\n$seq_right\n";
				}
			}
		}
	}
	close(TEMP);
	undef(%all2);
	undef(%all1);
	$self={
		range=>\%range
	};
	
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

